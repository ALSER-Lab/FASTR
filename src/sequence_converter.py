import os
import numpy as np 
import time 
import argparse
import cProfile
import pstats
import re
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from functools import partial
import io
import traceback


@dataclass
class FASTQRecord:
    index: int
    header: bytes
    sequence: np.ndarray
    quality: np.ndarray
    pair_number: int = 0

def detect_pair_number(header: str) -> int:
    """
    Detect if a read is pair 1 or pair 2 from its header (for paired-end compatability).
    Common patterns in paired-end fastq files are: /1, /2, or _1, _2
    """
    # Check for /1 or /2 pattern
    if '/1' in header:
        return 1
    elif '/2' in header:
        return 2
    elif '_1' in header:
        return 1
    elif '_2' in header:
        return 2
    
    return 0


ILLUMINA_PATTERN = re.compile(r'@([^:]+):(\d+):([^:]+):(\d+):(\d+):(\d+):(\d+)\s+(\d+):([YN]):(\d+):(.*)')
PACBIO_PATTERN = re.compile(r'@([^/]+)/(\d+)/(\d+)_(\d+)')
SRR_PATTERN = re.compile(r'@([A-Z]+)(\d+)\.(\d+)\s+(\d+)')

# Following functions extract IDs from respective machine headers AND common metadata: 
def parse_illumina_header(header: str) -> Tuple[Dict, str]:
    match = ILLUMINA_PATTERN.match(header)
    if not match:
        return {}, header  # Return as-is if doesn't match
    
    groups = match.groups()
    
    # Common metadata
    common = {
            'instrument': groups[0],
            'run_id': groups[1],
            'flowcell': groups[2],
            'lane': groups[3],
            }
    
    # Unique identifier: tile:x:y:read:filtered:control:barcode
    unique_id = f"{groups[4]}:{groups[5]}:{groups[6]}:{groups[7]}:{groups[8]}:{groups[9]}:{groups[10]}"
    return common, unique_id

def parse_pacbio_header(header: str) -> Tuple[Dict, str]:
    match = PACBIO_PATTERN.match(header)
    if not match:
        return {}, header
    
    common = {'movie': match.group(1)}
    # Unique identifier: zmw/start_end
    unique_id = f"{match.group(2)}/{match.group(3)}_{match.group(4)}"
    return common, unique_id

def parse_ont_header(header: str) -> Tuple[Dict, str]:
    parts = header.strip().split()
    if not parts:
        return {}, header
    
    read_id = parts[0][1:] if parts[0].startswith('@') else parts[0]  # Remove @
    
    # Parse key-value pairs
    kvs = {}
    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            kvs[key] = value
    
    # Common metadata
    common_keys = ['runid', 'sampleid', 'model_version_id', 'basecall_model_version_id']
    common = {k: kvs[k] for k in common_keys if k in kvs}
    
    # Extract unique fields
    unique_keys = ['read', 'ch', 'start_time']
    unique_parts = [f"{k}={kvs[k]}" for k in unique_keys if k in kvs]
    
    if unique_parts:
        unique_id = f"{read_id}:{':'.join(unique_parts)}"
    else:
        unique_id = read_id
    
    return common, unique_id


def parse_srr_header(header: str) -> Tuple[Dict, str]:
    match = SRR_PATTERN.match(header)
    if not match:
        return {}, header 
    
    prefix, accession, read_index, spot = match.groups()

    # Common metadata 
    common = {'prefix': prefix, 'accession': accession}
    unique_id = f"{read_index}:{spot}"

    return common, unique_id

def compress_header(header: str, sequencer_type: str) -> Tuple[Dict, str]:
    """
    Compress a single header based on sequencer type, returning the common metadata and unique id portion
    Returns (common_metadata_dict, unique_id)
    """
    if sequencer_type == 'illumina':
        return parse_illumina_header(header)
    elif sequencer_type == 'pacbio':
        return parse_pacbio_header(header)
    elif sequencer_type == 'ont':
        return parse_ont_header(header)
    elif sequencer_type == 'srr':
        return parse_srr_header(header)
    else:  # none or unrecognized
        return {}, header  # No compression

def format_metadata_header(common_metadata: Dict, sequencer_type: str) -> str:
    # Based on type of sequencer machine inputted returns formatted common metadata header
    if sequencer_type == 'illumina':
        return f"{common_metadata.get('instrument', '')}:{common_metadata.get('run_id', '')}:{common_metadata.get('flowcell', '')}:{common_metadata.get('lane', '')}"
    elif sequencer_type == 'pacbio':
        return common_metadata.get('movie', '')
    elif sequencer_type == 'ont':
        parts = []
        for key in ['runid', 'sampleid', 'model_version_id', 'basecall_model_version_id']:
            if key in common_metadata:
                parts.append(f"{key}={common_metadata[key]}")
        return ':'.join(parts)
    elif sequencer_type == 'srr':
        return f"{common_metadata.get('prefix', '')}{common_metadata.get('accession', '')}"
    else:
        return ''

def metadata_dict_equals(dict1: Dict, dict2: Dict) -> bool:
    """Compare two metadata dictionaries for equality, determining whether or not to flag a new flowcell in metadata"""
    if set(dict1.keys()) != set(dict2.keys()):
        return False
    for key in dict1:
        if dict1[key] != dict2[key]:
            return False
    return True

def get_scaling_equation(scaling_method: str, custom_formula=None, phred_alphabet_max=41) -> str:
    """
    Returns the mathematical equation string for each scaling method
    """
    if scaling_method == 'custom' and custom_formula:
        return custom_formula
    elif scaling_method == 'none':
        return "x"
    elif scaling_method == 'log':
        return f"1+62*(ln(x-1)/ln({phred_alphabet_max-1}))"
    elif scaling_method == 'linear':
        return "1+62*(x-40)/53"
    else:
        return "x"


def create_phred_quality_map(phred_offset=33, phred_alphabet_max=41): 
    # Offset is more of a custom thing for the user? It's not really in use, as the offset is applied during fastq gen itself
    phred_map = np.zeros(128, dtype=np.uint8)
    
    for ascii_val in range(128):
        quality_score = ascii_val - phred_offset
        # Clip to range [0, (phred_alphabet_max - 1)]
        quality_score = max(0, min(quality_score, (phred_alphabet_max-1)))
        phred_map[ascii_val] = quality_score
    
    return phred_map


def parse_custom_formula(formula: str, quality_scores: np.ndarray) -> np.ndarray:
    """
    Parse and evaluate a custom mathematical formula.
    Supports: +, -, *, /, **, (), ln(), log(), log10(), exp(), sqrt(), abs(), min(), max()
    Variable 'x' represents quality scores
    
    Example formulas are:
      "1 + 62 * (x - 40) / 53" - Linear
      "ln(x) * 10" - Weird
      "x ** 2 / 100" - Weird again
      "1 + 62 * (ln(x - 39) / ln(54))" - Log
    """
    # Remove f(x)= prefix if present
    formula = re.sub(r'^\s*f\s*\(\s*x\s*\)\s*=\s*', '', formula.strip())
    
    # Create safe namespace w/ numpy functions
    safe_dict = {'x': quality_scores, 'ln': np.log, 'log': np.log, 'log10': np.log10, 'exp': np.exp,
        'sqrt': np.sqrt, 'abs': np.abs, 'min': np.minimum, 'max': np.maximum, 'np': np, '__builtins__': {}}
    
    try:
        # Replace common math notation
        formula = formula.replace('^', '**')
        
        # Eval formula
        with np.errstate(divide='ignore', invalid='ignore'):
            result = eval(formula, safe_dict)
        
        # Handle scalar results (broadcast to array)
        if np.isscalar(result):
            result = np.full_like(quality_scores, result, dtype=np.float32)
        
        return result.astype(np.float32)
    
    except Exception as e:
        print(f"ERROR: Failed to parse custom formula '{formula}'")
        print(f"Error details: {e}")
        print("\nSupported operations: +, -, *, /, ** (power), (), ln(), log(), log10(), exp(), sqrt(), abs()")
        print("Variable 'x' represents quality scores")
        print("\nExample formulas:")
        print("  '1 + 62 * (x - 40) / 53'")
        print("  'ln(x - 39) / ln(54) * 62 + 1'")
        print("  'x ** 2 / 100'")
        raise

# Application of quality to bases
def apply_quality_to_bases(base_values, quality_scores, base_map, scaling_method='none', 
                           custom_formula=None, phred_alphabet_max=41):
    # Straight up one-hot encoding
    if scaling_method == 'none':
        return base_values
    
    # Calculate scale factors
    # Rather than individual calculations, we just calculate for all values then perform a dict lookup where needed 

    if scaling_method == 'custom':
        if custom_formula is None:
            print("ERROR: custom scaling method requires --custom_formula argument")
            return base_values
        
        # Parse and evaluate custom formula
        scale_factors = parse_custom_formula(custom_formula, quality_scores)

    elif scaling_method == 'log':
        # This makes larger quality scores closer in difference, and smaller quality scores more different from eachother. 
        # f(x) = 1 + 62 * (log(x-1) / log((phred_alphabet_max-1))), function works nicely by using a log10 curve (which the normal q score formula uses anyways)
        # Normal np implementation is: scale_factors = 1 + 62 * (np.log(quality_scores - 39) / np.log(54))
        # Do this to shut up divison by zero errors or other invalid vals that get clipped anyways
        with np.errstate(divide='ignore', invalid='ignore'): 
            LOG_DICT = np.clip(1 + 62 * (np.log(np.arange(0, 94, dtype=np.float32) - 1) / np.log((phred_alphabet_max-1))),1, 63).astype(np.uint8)
            scale_factors = LOG_DICT[quality_scores]

    elif scaling_method == 'linear':
        # Linear scaling
        scale_factors = 1 + 62 * (quality_scores - 40) / 53
    elif scaling_method == 'adaptive':
        # Linear scaling with min/max as bounds 
        min_q = quality_scores.min()
        max_q = quality_scores.max()
        if max_q == min_q:
            return np.full_like(base_values, 32, dtype=np.uint8)
        scale_factors = 1 + 62 * (quality_scores - min_q) / (max_q - min_q)
    else:
        return base_values
    
    # Clip and convert to uint8
    scale_factors = np.clip(scale_factors, 1, 63).astype(np.uint8)
    
    # Lookup table approach (fastest)
    base_min_lookup = np.zeros(256, dtype=np.uint8)
    base_min_lookup[63] = 1      # A
    base_min_lookup[127] = 65    # T
    base_min_lookup[191] = 129   # C
    base_min_lookup[255] = 193   # G
    base_min_lookup[1] = 0       # N, color of 1 though I should probably change this 
    
    result = base_min_lookup[base_values] + scale_factors - 1
    
    return result.astype(np.uint8)

def process_chunk_worker(chunk_data, base_map, phred_map, compress_headers, 
                        sequencer_type, paired_end, keep_bases, binary_bases, 
                        keep_quality, binary_quality, quality_scaling, custom_formula,
                        phred_alphabet_max, min_quality, BYTE_LOOKUP, binary,
                        remove_repeating_header):
    """
    Worker function method that processes a single chunk in parallel.
    """
    try:
        chunk_id, buffer, start_index = chunk_data
        print(f"Worker processing chunk {chunk_id} with {len(buffer)} bytes")  # debug 
        
        # parse the records
        records, _, metadata, count = parse_fastq_records_from_buffer(
            buffer, start_index, base_map, phred_map,
            compress_headers, sequencer_type, paired_end,
            keep_bases, binary_bases, keep_quality, binary_quality
        )
        
        print(f"Worker parsed {count} records from chunk {chunk_id}")  # debug
        
        # Process records and write to in mem buffer
        output_buffer = io.BytesIO()
        
        if records:
            process_and_write_records(
                records, output_buffer, base_map, quality_scaling,
                custom_formula, phred_alphabet_max, min_quality, keep_bases,
                binary_bases, binary, keep_quality,
                remove_repeating_header, compress_headers, BYTE_LOOKUP
            )
        
        print(f"Worker completed chunk {chunk_id}")  # debug
        return (chunk_id, output_buffer.getvalue(), metadata, count)
    
    except Exception as e: # I don't know HOW an error would arise honestly, but this is here just for safety
        print(f"ERROR in worker processing chunk {chunk_id}: {e}")
        traceback.print_exc()
        raise


def parse_fastq_records_from_buffer(buffer: bytes, start_index: int, base_map: np.ndarray, 
                                    phred_map: Optional[np.ndarray], compress_headers: bool, 
                                    sequencer_type: str, paired_end: bool, keep_bases: bool, 
                                    binary_bases: bool, keep_quality: bool, binary_quality: bool) -> Tuple[List[FASTQRecord], bytes, Optional[Dict], int]:
    
    records = []
    current_flowcell_metadata = None
    lines = buffer.split(b'\n') # Split entire buffer ONCE
    num_complete_records = (len(lines) - 1) // 4 # Calculate amt of complete FASTQ records
    
    if num_complete_records == 0: # Quick return
        return records, buffer, current_flowcell_metadata, 0
    
    
    leftover = b'\n'.join(lines[num_complete_records * 4:]) # Get amount before we process
    
    # Pre-create binary map ONCE if needed
    binary_map = None
    if binary_bases:
        binary_map = np.zeros(128, dtype=np.uint8)
        binary_map[ord('A')] = 0
        binary_map[ord('T')] = 1
        binary_map[ord('C')] = 2
        binary_map[ord('G')] = 3
        binary_map[ord('N')] = 4
    
    # Process in smaller batches to avoid numpy split overhead BUT big enough to amortize loop costs
    BATCH_SIZE = 10000  # Process 10k records at a time (this val honestly doesn't even really matter too much)
    
    for batch_start in range(0, num_complete_records, BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, num_complete_records)
        batch_size = batch_end - batch_start
        
        line_start = batch_start * 4
        line_end = batch_end * 4
        
        
        batch_lines = lines[line_start:line_end] # Get batch lines
        # Vectorized operations on this batch! Ilysm numpy
        header_strs = [batch_lines[i].decode('utf-8', errors='ignore') for i in range(0, len(batch_lines), 4)]
        seq_data_list = [batch_lines[i].replace(b'\r', b'') for i in range(1, len(batch_lines), 4)]
        qual_data_list = [batch_lines[i].replace(b'\r', b'') for i in range(3, len(batch_lines), 4)]
        
        seq_lengths = [len(s) for s in seq_data_list] # Concatenate and convert ONCE for entire batch (major speedup)
        all_seq_bytes = b''.join(seq_data_list)
        all_seq_array = np.frombuffer(all_seq_bytes, dtype=np.uint8) # Single frombuffer for entire batch (major speedup)
        
        # Apply base mapping to ENTIRE batch at once (avoid overhead from per-record way)
        if keep_bases:
            mapped_batch = all_seq_array
        elif binary_map is not None:
            mapped_batch = binary_map[all_seq_array]
        else:
            mapped_batch = base_map[all_seq_array]
        
        cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
        
        # Create views (instead of copies, which led to a major speedup as well)
        seq_arrays = [mapped_batch[cumsum_lengths[i]:cumsum_lengths[i+1]] for i in range(batch_size)]
        
        
        qual_arrays = [] # Process qualities (if needed)
        if phred_map is not None:
            all_qual_bytes = b''.join(qual_data_list)
            all_qual_array = np.frombuffer(all_qual_bytes, dtype=np.uint8)
            qual_mapped = phred_map[all_qual_array]
            # Use views for quality too
            qual_arrays = [qual_mapped[cumsum_lengths[i]:cumsum_lengths[i+1]] for i in range(batch_size)]
        
        # Build records
        for idx in range(batch_size):
            record_index = start_index + batch_start + idx
            header_str = header_strs[idx]
            
            # A fast path is to skip pair detection if not needed
            pair_number = 0
            if paired_end:
                # Inline detection (avoids function call, to squeeze some speedups)
                if '/1' in header_str or '_1' in header_str:
                    pair_number = 1
                elif '/2' in header_str or '_2' in header_str:
                    pair_number = 2
            
            # Header compression
            if compress_headers and sequencer_type != 'none':
                common, unique_id = compress_header(header_str, sequencer_type)
                if common:
                    current_flowcell_metadata = common
                
                if paired_end and pair_number > 0:
                    header = f"@{unique_id}/{pair_number}\n".encode('utf-8')
                else:
                    header = f"@{unique_id}\n".encode('utf-8')
            else:
                header = batch_lines[idx * 4] + b'\n'
            
            # Get quality
            qual_vals = qual_arrays[idx] if qual_arrays else np.array([], dtype=np.uint8)
            
            record = FASTQRecord(
                index=record_index,
                header=header,
                sequence=seq_arrays[idx],
                quality=qual_vals,
                pair_number=pair_number
            )
            records.append(record)
    
    return records, leftover, current_flowcell_metadata, num_complete_records

def process_and_write_records(records: List[FASTQRecord], outfile, base_map: np.ndarray,
                               quality_scaling: str, custom_formula: Optional[str],
                               phred_alphabet_max: int, min_quality: int, keep_bases: bool,
                               binary_bases: bool, binary: bool, keep_quality: bool,
                               remove_repeating_header: bool, compress_headers: bool,
                               BYTE_LOOKUP: np.ndarray):
    """
    This function processes quality scaling on records and then immediately writes to file.
    """
    if not records:
        return
    
    seq_lengths = [len(r.sequence) for r in records]  # Store sequence lengths before concat
    flat_mapped = np.concatenate([r.sequence for r in records]).astype(np.uint8)
    
    all_quals = [r.quality for r in records if len(r.quality) > 0]
    
    if len(all_quals) > 0 and quality_scaling != 'none':
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        
        if not keep_bases and not binary_bases:
            # Apply scaling on flattened arrays 
            flat_mapped = apply_quality_to_bases(
                flat_mapped, flat_quals, base_map, quality_scaling,
                custom_formula, phred_alphabet_max
            )
    
    if min_quality > 0 and len(all_quals) > 0:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        low_quality_mask = flat_quals < min_quality
        
        if keep_bases:
            flat_mapped[low_quality_mask] = ord('N')  # Keep as ASCII
        elif binary_bases:
            flat_mapped[low_quality_mask] = 4  # N = 4 in binary encoding
        else:
            flat_mapped[low_quality_mask] = base_map[ord('N')]  # Convert to numeric

    cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
    processed_seqs = [
        flat_mapped[cumsum_lengths[i]:cumsum_lengths[i+1]]
        for i in range(len(records))
    ]
    
    # Write all sequences
    for idx, record in enumerate(records):
        # Only write headers if remove_repeating_header is disabled
        if not remove_repeating_header:
            if compress_headers:
                outfile.write(record.header)
            else:
                outfile.write(record.header)
        
        outfile.write(b'<')  # Add start marker for sequence
        
        seq_data = processed_seqs[idx]
        
        if binary:
            outfile.write(seq_data.tobytes())
        else:
            if keep_bases:
                outfile.write(seq_data.tobytes())  # Write ASCII directly
            elif binary_bases:
                outfile.write(b''.join(BYTE_LOOKUP[seq_data]))  # Write binary bases as text numbers for inspection
            else:
                outfile.write(b''.join(BYTE_LOOKUP[seq_data]))  # Convert numbers to text
        
        outfile.write(b'\n')
        
        # Write quality scores if keep_quality is enabled
        if keep_quality and len(record.quality) > 0:
            outfile.write(b'+\n')
            outfile.write(record.quality)  # ASCII quality string
            outfile.write(b'\n')



def export_scalars_to_txt(fastq_path, base_map, output_path, phred_map=None, min_quality=0, 
                          quality_scaling='none', binary=True, log_a=None, compress_headers=False, 
                          sequencer_type='none', sra_accession=None, keep_bases=False, keep_quality=False, 
                          binary_bases=False, binary_quality=False, custom_formula=None, multiple_flowcells=False,
                          remove_repeating_header=False, phred_alphabet_max=41, paired_end=False, 
                          paired_end_mode='same_file', chunk_size_mb=32, num_workers=4):
    """
    This function uses a streaming architecture to handle files of any size w/ constant mem usage.
    This is where we write the binary vals, and this function acts as a "main" for our base scaling, parsing etc.
    Uses streaming architecture with pool.imap (instead of normal pool.map) to try and maintain low (not low, moreso consistent) memory usage.
    """
    
    file_size = os.path.getsize(fastq_path)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    print(f"Using {num_workers} worker processes for parallel processing")
    
    if compress_headers and sequencer_type != 'none':
        print(f"Header compression enabled (sequencer type: {sequencer_type})")
        if multiple_flowcells:
            print("Multiple flowcell detection enabled")
        if remove_repeating_header:
            print("Repeating header metadata removal enabled - only unique IDs will be stored")
    
    if paired_end:
        print(f"Paired-end mode enabled (output mode: {paired_end_mode})")
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    # Premake lookup table
    BYTE_LOOKUP = np.array([f'{i:03d}'.encode('ascii') for i in range(256)])
    
    # For tracking flowcells and sequences
    flowcell_metadata_list = []
    current_flowcell_metadata = None
    flowcell_start_index = 0
    total_sequences = 0
    
    print("Reading and processing chunks in parallel...")
    
    def chunk_generator(): # Generator function to yield chunks as they're read
        buffer = b''
        chunk_id = 0
        start_index = 0
        
        with open(fastq_path, 'rb', buffering=chunk_size_bytes) as infile:
            while True:
                chunk = infile.read(chunk_size_bytes)
                if not chunk and not buffer:
                    break
                
                buffer += chunk
                
                last_at = buffer.rfind(b'\n@') # Find last complete record boundary
                if last_at == -1 or not chunk:
                    process_buffer = buffer
                    buffer = b''
                else:
                    process_buffer = buffer[:last_at+1]
                    buffer = buffer[last_at+1:]
                
                if process_buffer:
                    yield (chunk_id, process_buffer, start_index) # Yield chunk data for processing
                    
                    # estmate number of records for next chunk's start index
                    estimated_records = process_buffer.count(b'\n@')
                    start_index += estimated_records
                    chunk_id += 1
                
                if chunk_id % 10 == 0: # maybe this isn't wise right now? I'll tweak this on deployment to not barf out a ton of text
                    print(f"Read {chunk_id} chunks...")
    
    with open(output_path, 'wb', buffering=chunk_size_bytes) as outfile:
        # Write SRA accession if provided
        if sra_accession:
            outfile.write(f"#{sra_accession}\n".encode('utf-8'))
        
        # Reserve space for metadata headers (we'll come back to write these)
        metadata_position = outfile.tell()
        MAX_METADATA_LINES = 20  
        placeholder_size = MAX_METADATA_LINES * 200  # 200 chars per line should be enough
        if compress_headers and sequencer_type != 'none':
            outfile.write(b' ' * placeholder_size)  # Write spaces as placeholder
            outfile.write(b'\n')
        
        # create worker pool and process chunks as they come in
        with Pool(processes=num_workers) as pool:
            # Create partial function w/ fixed args
            worker_func = partial(
                process_chunk_worker,
                base_map=base_map,
                phred_map=phred_map,
                compress_headers=compress_headers,
                sequencer_type=sequencer_type,
                paired_end=paired_end,
                keep_bases=keep_bases,
                binary_bases=binary_bases,
                keep_quality=keep_quality,
                binary_quality=binary_quality,
                quality_scaling=quality_scaling,
                custom_formula=custom_formula,
                phred_alphabet_max=phred_alphabet_max,
                min_quality=min_quality,
                BYTE_LOOKUP=BYTE_LOOKUP,
                binary=binary,
                remove_repeating_header=remove_repeating_header
            )
            
            # Use imap to process chunks as they're read (streaming) -- previous implementation used normal pool.map, which loaded the file into memory
            # chunksize=1 makes sure order is preserved
            for chunk_id, processed_bytes, metadata, count in pool.imap(worker_func, chunk_generator(), chunksize=1):
                # Write immediately as each chunk completes
                outfile.write(processed_bytes)
                
                if metadata: # Track flowcell metadata
                    if multiple_flowcells:
                        if current_flowcell_metadata is None:
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                        elif not metadata_dict_equals(current_flowcell_metadata, metadata):
                            flowcell_metadata_list.append((
                                current_flowcell_metadata.copy(),
                                flowcell_start_index,
                                total_sequences - 1
                            ))
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                            print(f"Flowcell change detected at sequence {total_sequences}")
                    elif current_flowcell_metadata is None:
                        current_flowcell_metadata = metadata
                
                total_sequences += count
                
                if total_sequences % 100000 == 0:
                    print(f"Processed {total_sequences:,} sequences...")
        
        # Finalize flowcell tracking
        if multiple_flowcells and current_flowcell_metadata is not None:
            flowcell_metadata_list.append((
                current_flowcell_metadata.copy(),
                flowcell_start_index,
                total_sequences - 1
            ))
        elif current_flowcell_metadata:
            flowcell_metadata_list.append((current_flowcell_metadata, 0, total_sequences - 1))
        
        # Print flowcell summary (when we start testing it later on)
        if multiple_flowcells and len(flowcell_metadata_list) > 1:
            print(f"\nDetected {len(flowcell_metadata_list)} flowcells:")
            for idx, (metadata, start, end) in enumerate(flowcell_metadata_list):
                fc_id = metadata.get('flowcell', metadata.get('movie', 'unknown'))
                print(f"  Flowcell {idx + 1}: {fc_id} (sequences {start}-{end}, total: {end - start + 1})")
        
        # Now go back and write metadata headers at the beginning 
        # This is why when you try to open the file while data is being written, the first 1-2 line(s) are blank! We write common metadata after everything else
        if compress_headers and sequencer_type != 'none':
            end_position = outfile.tell()
            outfile.seek(metadata_position)
            
            metadata_lines = []
            
            if multiple_flowcells and len(flowcell_metadata_list) > 0:
                # Build metadata lines
                for metadata, start_idx, end_idx in flowcell_metadata_list:
                    metadata_line = format_metadata_header(metadata, sequencer_type)
                    metadata_lines.append(f"#{metadata_line}\n")
                    metadata_lines.append(f"#SEQUENCER:{sequencer_type}\n")
                    equation = get_scaling_equation(quality_scaling, custom_formula, phred_alphabet_max)
                    metadata_lines.append(f"#{equation}\n")
                    metadata_lines.append(f"#RANGE:{start_idx}-{end_idx}\n")
            else:
                # Single flowcell
                if current_flowcell_metadata:
                    metadata_line = format_metadata_header(current_flowcell_metadata, sequencer_type)
                    metadata_lines.append(f"#{metadata_line}\n")
                metadata_lines.append(f"#SEQUENCER:{sequencer_type}\n")
                equation = get_scaling_equation(quality_scaling, custom_formula, phred_alphabet_max)
                metadata_lines.append(f"#{equation}\n")
            
            # Write metadata and pad remaining space
            metadata_bytes = ''.join(metadata_lines).encode('utf-8')
            outfile.write(metadata_bytes)
            remaining = placeholder_size - len(metadata_bytes)
            if remaining > 0:
                outfile.write(b' ' * remaining)
            outfile.write(b'\n')
            
            outfile.seek(end_position)  # Return to end
    
    print(f"Total sequences written: {total_sequences:,}")

def main():
    argument_parser = argparse.ArgumentParser(description="Convert and compress FASTQ/FASTA files to scalar format")
    argument_parser.add_argument("input_path", type=str, help="Path of .fasta or .fastq file")
    argument_parser.add_argument("output_path", type=str, help="Output file path")

    # Quick mode implementation
    argument_parser.add_argument("--mode", type=int,
                            help="What mode to use when writing out converted values." \
                            "1: Header compression only" \
                            "2: Base conversion into numbers only" \
                            "3: Header and base conversion, written out in two lines" \
                            "4: Repeating header removal entirely, base conversion kept, written out in one line")
    # Paired-end args
    argument_parser.add_argument("--paired_end", type=int, default=0,
                                help="Signifies if the input file contains paired-end reads (0/1, default 0)")
    argument_parser.add_argument("--paired_end_mode", type=str, 
                                 choices=['same_file', 'separate_files'],
                                 default='same_file',
                                 help="Output mode for paired-end reads: same_file or separate_files (default: same_file)")
    
    # SRA accession
    argument_parser.add_argument("--sra_accession", type=str, default=None,
                                help="SRA accession number (e.g., SRR12345678)")
    
    # Header compression arguments
    argument_parser.add_argument("--compress_headers", type=int, default=0,
                                help="Compress FASTQ headers on-the-fly (0/1, default 0)")
    argument_parser.add_argument("--sequencer_type", type=str, 
                                choices=['none', 'illumina', 'pacbio', 'ont', 'srr'],
                                default='none',
                                help="Sequencer type for header compression (default: none)")
    argument_parser.add_argument("--multiple_flowcells", type=int, default=0,
                                help="Enable multiple flowcell detection and tracking (0/1, default 0)")
    
    # Quality arguments
    argument_parser.add_argument("--extract_quality", type=int, default=1, 
                                help="For FASTQ: extract quality scores (0/1, default 1)")
    argument_parser.add_argument("--phred_offset", type=int, default=33,
                                help="Phred quality offset (if needed)(default 33)")
    argument_parser.add_argument("--min_quality", type=int, default=0,
                                help="Minimum quality score threshold (default 0)")
    argument_parser.add_argument("--quality_scaling", type=str, 
                                 choices=['log','log_custom', 'custom'], 
                                 default='none',
                                 help="Quality scaling method (default: none)")
    argument_parser.add_argument("--custom_formula", type=str, default=None,
                                help="Custom formula for quality scaling (use 'x' for quality score). "
                                     "Example: '1 + 62 * (x - 40) / 53' or 'ln(x) * 10'")
    argument_parser.add_argument("--log_a", type=int, default=0,
                                help="Tunable a value for application in custom log function (default 0)")
    # Output format
    argument_parser.add_argument("--binary_write", type=int, default=1,
                                help="Enable binary writing of sequence integers (0/1, default 1)")
    argument_parser.add_argument("--keep_bases", type=int, default=0,
                                help="Whether or not to return textual bases w/o any scaling or one-hot encoding (0/1, default 0)")
    argument_parser.add_argument("--keep_quality", type=int, default=0,
                            help="Keep original quality scores in output (0/1, default 0)")
    argument_parser.add_argument("--binary_bases", type=int, default=0,
                            help="Use binary encoding for STRING BASES, for binary writing of scaled integers use --binary_write 1 (0/1, default 0)")
    argument_parser.add_argument("--binary_quality", type=int, default=0,
                            help="Write quality scores as binary numeric values instead of ASCII (0/1, default 0)")
    
    argument_parser.add_argument("--remove_repeating_header", type=int, default=0,
                            help="Remove repeating metadata from individual sequence headers, store only at top (0/1, default 0)")
    argument_parser.add_argument("--phred_alphabet", type=str, default='phred42',
                            help="What phred quality (q-score) ascii character alphabet is used by the inputted fastq")

    # Base mapping
    argument_parser.add_argument("--gray_N", type=int, default=1, help="Grayscale value for N (default 1)")
    argument_parser.add_argument("--gray_A", type=int, default=63, help="Grayscale value for A (default 63)")
    argument_parser.add_argument("--gray_T", type=int, default=127, help="Grayscale value for T (default 127)")
    argument_parser.add_argument("--gray_C", type=int, default=191, help="Grayscale value for C (default 191)")
    argument_parser.add_argument("--gray_G", type=int, default=255, help="Grayscale value for G (default 255)")
    
    # Profiling
    argument_parser.add_argument("--profile", type=int, default=0,
                                help="Enable profiling (0/1, default 0)")
    
    # Multiprocessing
    argument_parser.add_argument("--num_workers", type=int, default=1,
                                help="Number of parallel workers (default: 1, use 4+ for large files >5GB)")

    user_arguments = argument_parser.parse_args()
    alphabet_max = None

    if user_arguments.keep_bases == 1 and user_arguments.binary_bases == 1:
        print("ERROR: Cannot use both --keep_bases and --binary_bases. Choose one:")
        print("  --keep_bases 1: Keep ASCII letters (A=65, T=84, etc.)")
        print("  --binary_bases 1: Convert to compact encoding (A=0, T=1, etc.)")
        exit(1)
    
    if user_arguments.quality_scaling == 'custom' and user_arguments.custom_formula is None:
        print("ERROR: --quality_scaling custom requires --custom_formula argument")
        print("\nExample usage:")
        print("  --quality_scaling custom --custom_formula '1 + 62 * (x - 40) / 53'")
        print("  --quality_scaling custom --custom_formula 'ln(x - 39) / ln(54) * 62 + 1'")
        print("  --quality_scaling custom --custom_formula 'sqrt(x) * 5'")
        exit(1)
    
    if user_arguments.multiple_flowcells == 1 and user_arguments.compress_headers == 0:
        print("WARNING: --multiple_flowcells requires --compress_headers 1 to function")
        print("Enabling header compression automatically...")
        user_arguments.compress_headers = 1
    
    if user_arguments.remove_repeating_header == 1 and user_arguments.compress_headers == 0:
        print("WARNING: --remove_repeating_header requires --compress_headers 1 to function")
        print("Enabling header compression automatically...")
        user_arguments.compress_headers = 1

    if user_arguments.mode == 1:
        user_arguments.compress_headers = 1
        user_arguments.binary_write = 0

    elif user_arguments.mode == 2:
        user_arguments.compress_headers = 0
        user_arguments.binary_write = 1

    elif user_arguments.mode == 3:
        user_arguments.compress_headers = 1
        user_arguments.binary_write = 1
    
    elif user_arguments.mode == 4:
        user_arguments.compress_headers = 1
        user_arguments.remove_repeating_header = 1
        user_arguments.binary_write = 1

    if user_arguments.phred_alphabet == "phred42":
        phred_alphabet_max = 41
    elif user_arguments.phred_alphabet == "phred63":
        phred_alphabet_max = 62
    elif user_arguments.phred_alphabet == "phred94":
        phred_alphabet_max = 93
    

    # Create base map
    numpy_base_map = np.zeros(128, dtype=np.uint8) 
    numpy_base_map[ord("N")] = user_arguments.gray_N
    numpy_base_map[ord("A")] = user_arguments.gray_A
    numpy_base_map[ord("T")] = user_arguments.gray_T
    numpy_base_map[ord("C")] = user_arguments.gray_C
    numpy_base_map[ord("G")] = user_arguments.gray_G

    start_time = time.perf_counter()

    # Initialize profiler if requested (for dev optim)
    profiler = None
    if user_arguments.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        print("Profiling enabled...")

    # Create phred map if needed
    phred_map = create_phred_quality_map(user_arguments.phred_offset, phred_alphabet_max) if user_arguments.extract_quality else None

    print(f"Converting sequences from {user_arguments.input_path}...")
    
    export_scalars_to_txt(
        user_arguments.input_path, numpy_base_map, 
        user_arguments.output_path, phred_map, user_arguments.min_quality, 
        user_arguments.quality_scaling, user_arguments.binary_write,
        user_arguments.log_a, compress_headers=(user_arguments.compress_headers == 1),
        sequencer_type=user_arguments.sequencer_type,
        sra_accession=user_arguments.sra_accession,
        keep_bases=user_arguments.keep_bases, keep_quality=user_arguments.keep_quality,
        binary_bases=user_arguments.binary_bases, binary_quality=user_arguments.binary_quality,
        custom_formula=user_arguments.custom_formula,
        multiple_flowcells=(user_arguments.multiple_flowcells == 1),
        remove_repeating_header=(user_arguments.remove_repeating_header == 1),
        phred_alphabet_max=phred_alphabet_max,
        paired_end=(user_arguments.paired_end == 1),
        paired_end_mode=user_arguments.paired_end_mode,
        num_workers=user_arguments.num_workers
        )

    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Task completed in {elapsed_time:.4f} seconds")
    print(f"Saved scalars to {user_arguments.output_path}")
    
    if profiler is not None:
        profiler.disable()
        print("\n" + "="*80)
        print("Profiling Results:")
        print("="*80)
        stats = pstats.Stats(profiler)
        stats.sort_stats('cumulative')
        stats.print_stats(30)
        stats.sort_stats('tottime')
        print("\n" + "="*80)
        print("By total time:")
        print("="*80)
        stats.print_stats(30)
if __name__ == "__main__":
    main()