import os
import numpy as np 
import time 
import argparse
import cProfile
import pstats
import re
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

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

# Following functions extract IDs from respective machine headers AND common metadata: 
def parse_illumina_header(header: str) -> Tuple[Dict, str]:
    pattern = re.compile(r'@([^:]+):(\d+):([^:]+):(\d+):(\d+):(\d+):(\d+)\s+(\d+):([YN]):(\d+):(.*)')
    match = pattern.match(header)
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
    pattern = re.compile(r'@([^/]+)/(\d+)/(\d+)_(\d+)')
    match = pattern.match(header)
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

def get_scaling_equation(scaling_method: str,min_q=None, max_q=None, custom_formula=None, phred_alphabet_max=41) -> str:
    """
    Returns the mathematical equation string for each scaling method
    """
    if scaling_method == 'custom' and custom_formula:
        return custom_formula
    elif scaling_method == 'none':
        return "x"
    elif scaling_method == 'log':
        return f"1+62*(ln(x-1)/ln({phred_alphabet_max-1}))"
    elif scaling_method == 'log_adaptive':
        min_val = min_q if min_q is not None else 40
        max_val = max_q if max_q is not None else 93
        return f"1+62*(ln(x-{min_val-1})/ln({max_val-min_val-1}))"
    elif scaling_method == 'linear':
        return "1+62*(x-40)/53"
    else:
        return "x"


def create_phred_quality_map(phred_offset=0, phred_alphabet_max=41): 
    # Offset is more of a custom thing for the user? It's not really in use, as the offset is applied during fastq gen itself
    phred_map = np.zeros(128, dtype=np.uint8)
    
    for ascii_val in range(128):
        quality_score = ascii_val - phred_offset
        # Clip to range [0, (phred_alphabet_max - 1)]
        quality_score = max(40, min(quality_score, (phred_alphabet_max-1)))
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
                           dataset_min_q=None, dataset_max_q=None, custom_formula=None, 
                           phred_alphabet_max=41):
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

    elif scaling_method == 'log_adaptive':
        # Use dataset-wide min/max passed in
        if dataset_min_q is None or dataset_max_q is None:
            print("ERROR: log_adaptive requires dataset min/max quality scores")
            return base_values
            
        min_q = dataset_min_q
        max_q = dataset_max_q
    
        # Avoid divide by zero if all values are the same
        if max_q == min_q:
            scale_factors = np.full_like(quality_scores, 32, dtype=np.uint8)
        else:
            x = np.arange(min_q, max_q + 1, dtype=np.float32)  # Only our relevant range
            # Tunable log curve mapping (min_q, 1) and (max_q, 63)
            LOG_ADP_DICT = 1 + 62 * (np.log(x - (min_q - 1)) / np.log(max_q - (min_q - 1)))
            LOG_ADP_DICT = np.clip(LOG_ADP_DICT, 1, 63).astype(np.uint8)

            scale_table = np.zeros(94, dtype=np.uint8)
            scale_table[min_q:max_q+1] = LOG_ADP_DICT
            scale_factors = scale_table[quality_scores]
 
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


def parse_fastq_records_from_buffer(buffer: bytes, start_index: int, base_map: np.ndarray, 
                                    phred_map: Optional[np.ndarray], compress_headers: bool, 
                                    sequencer_type: str, paired_end: bool, keep_bases: bool, 
                                    binary_bases: bool, keep_quality: bool, binary_quality: bool) -> Tuple[List[FASTQRecord], bytes, Optional[Dict], int]:
    """
    This method parses FASTQ records from a buffer and returns complete records/leftover incomplete data.
    """
    records = []
    current_flowcell_metadata = None
    sequence_count = 0
    i = 0
    
    while i < len(buffer):
        if buffer[i:i+1] == b'@':
            header_end = buffer.find(b'\n', i)
            if header_end == -1:
                # This means we have an incomplete record, so we'll just return what we have so far
                return records, buffer[i:], current_flowcell_metadata, sequence_count
            
            header_str = buffer[i:header_end].decode('utf-8', errors='ignore')
            
            # Detect pair num if paired-end arg is enabled
            pair_number = 0
            if paired_end:
                pair_number = detect_pair_number(header_str)
            
            # Compress header if enabled and extract metadata from first header
            header = None
            if compress_headers and sequencer_type != 'none':
                common, unique_id = compress_header(header_str, sequencer_type)
                if common:
                    current_flowcell_metadata = common
                
                # Append pair number to header for modes 1/3
                if paired_end and pair_number > 0:
                    header = f"@{start_index + sequence_count}:{unique_id}/{pair_number}\n".encode('utf-8')
                else:
                    header = f"@{start_index + sequence_count}:{unique_id}\n".encode('utf-8')
            else:
                header = buffer[i:header_end+1]
            # find base sequence
            seq_end = buffer.find(b'\n', header_end + 1)
            if seq_end == -1:
                return records, buffer[i:], current_flowcell_metadata, sequence_count
            
            seq_data = buffer[header_end+1:seq_end].replace(b'\r', b'')
            
            # find '+' line
            plus_end = buffer.find(b'\n', seq_end + 1)
            if plus_end == -1:
                return records, buffer[i:], current_flowcell_metadata, sequence_count
            
            # find quality line(s)
            qual_end = buffer.find(b'\n', plus_end + 1)
            if qual_end == -1:
                # check if we might have a complete record at end of buffer
                potential_qual_length = len(seq_data)
                if plus_end + 1 + potential_qual_length < len(buffer):
                    qual_end = plus_end + 1 + potential_qual_length
                else:
                    return records, buffer[i:], current_flowcell_metadata, sequence_count
            
            qual_data = buffer[plus_end+1:qual_end].replace(b'\r', b'')
            
            # conditional to check for corrupted data length (mismatch between base num and quality size)
            if len(seq_data) != len(qual_data):
                print(f"Skipping corrupted record {start_index + sequence_count} due to length mismatch: {len(seq_data)} bases vs {len(qual_data)} quality scores.")
                i = qual_end if qual_end < len(buffer) else len(buffer)
                continue # skip current record
            
            # Process quality WITHOUT min/max calculation
            qual_vals = None
            qual_string = None
            
            if phred_map is not None:
                qual_ascii = np.frombuffer(qual_data, dtype=np.uint8)
                qual_vals = phred_map[qual_ascii]
                
                # Store original quality string if keep_quality is enabled
                if keep_quality:
                    if binary_quality:
                        # Store the same numeric values we just calculated
                        qual_string = qual_vals.copy()  # Use .copy() to avoid reference issues
                    else:
                        qual_string = qual_data
            
            # Process sequence
            ascii_vals = np.frombuffer(seq_data, dtype=np.uint8)
            
            if keep_bases:
                mapped_vals = ascii_vals  # Keep original ASCII values 
            elif binary_bases:
                binary_map = np.zeros(128, dtype=np.uint8)
                binary_map[ord('A')] = 0
                binary_map[ord('T')] = 1
                binary_map[ord('C')] = 2
                binary_map[ord('G')] = 3
                binary_map[ord('N')] = 4
                mapped_vals = binary_map[ascii_vals]
            else:
                mapped_vals = base_map[ascii_vals]  # Convert to 255 range
            
            # store new fastq record
            record = FASTQRecord(
                index=start_index + sequence_count,
                header=header,
                sequence=mapped_vals,
                quality=qual_vals if qual_vals is not None else np.array([], dtype=np.uint8),
                pair_number=pair_number
            )
            records.append(record)
            sequence_count += 1
            
            i = qual_end
        else:
            i += 1
    
    return records, b'', current_flowcell_metadata, sequence_count

def process_and_write_records(records: List[FASTQRecord], outfile, base_map: np.ndarray,
                               quality_scaling: str, dataset_min_q: Optional[int], 
                               dataset_max_q: Optional[int], custom_formula: Optional[str],
                               phred_alphabet_max: int, min_quality: int, keep_bases: bool,
                               binary_bases: bool, binary: bool, keep_quality: bool,
                               remove_repeating_header: bool, BYTE_LOOKUP: np.ndarray):
    """
    This function processes quality scaling on records and then immediately writes to file.
    """
    if not records:
        return
    
    seq_lengths = [len(r.sequence) for r in records] # Store sequence lengths before concat
    flat_mapped = np.concatenate([r.sequence for r in records]).astype(np.uint8)
    
    all_quals = [r.quality for r in records if len(r.quality) > 0]
    
    if len(all_quals) > 0 and quality_scaling != 'none':
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        
        if not keep_bases and not binary_bases:
            # Apply scaling on flattened arrays 
            flat_mapped = apply_quality_to_bases(flat_mapped, flat_quals, base_map, quality_scaling, 
                                                 dataset_min_q, dataset_max_q, custom_formula, phred_alphabet_max)
    
    if min_quality > 0 and len(all_quals) > 0:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        low_quality_mask = flat_quals < min_quality
        
        if keep_bases:
            flat_mapped[low_quality_mask] = ord('N')  # Keep as ASCII
        elif binary_bases:
            flat_mapped[low_quality_mask] = 4  # N = 4 in binary encoding
        else:
            flat_mapped[low_quality_mask] = base_map[ord('N')]  # Convert to numeric

    split_indices = np.cumsum(seq_lengths)[:-1] # Split back into individual sequences
    processed_seqs = list(np.split(flat_mapped, split_indices))
    
    # Write all sequences
    for idx, record in enumerate(records):
        # Store only unique ID if remove_repeating_header is enabled
        if not remove_repeating_header:
            outfile.write(record.header)
        
        outfile.write(b'<')  # Add start marker for sequence (because sequences may be written out on multiple lines, sort of like FASTA's >'' signifier)
        if binary:
            outfile.write(processed_seqs[idx].tobytes())
        else:
            if keep_bases:
                outfile.write(processed_seqs[idx].tobytes())  # Write ASCII directly
            elif binary_bases:
                outfile.write(b''.join(BYTE_LOOKUP[processed_seqs[idx]])) # Write binary bases as text numbers for inspection
            else:
                outfile.write(b''.join(BYTE_LOOKUP[processed_seqs[idx]]))  # Convert numbers to text
        
        outfile.write(b'\n')
        
        # Write quality scores if keep_quality is enabled
        if keep_quality and len(record.quality) > 0:
            outfile.write(b'+\n')
            if binary_quality:
                outfile.write(record.quality.tobytes()) # Write as binary bytes (numeric quality values)
            else:
                # Write as ASCII text (original quality string)
                outfile.write(record.quality)
            outfile.write(b'\n')

def export_scalars_to_txt(fastq_path, base_map, output_path, phred_map=None, min_quality=0, 
                          quality_scaling='none', binary=True, log_a=None, compress_headers=False, 
                          sequencer_type='none', sra_accession=None, keep_bases=False, keep_quality=False, 
                          binary_bases=False, binary_quality=False, custom_formula=None, multiple_flowcells=False,
                          remove_repeating_header=False, phred_alphabet_max=41, paired_end=False, 
                          paired_end_mode='same_file', chunk_size_mb=100):
    """
    This function uses a streaming architecture to handle files of any size w/ constant mem usage.
    This is where we write the binary vals, and this function acts as a "main" for our base scaling, parsing etc. 
    """
    
    file_size = os.path.getsize(fastq_path)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    
    if compress_headers and sequencer_type != 'none':
        print(f"Header compression enabled (sequencer type: {sequencer_type})")
        if multiple_flowcells:
            print("Multiple flowcell detection enabled")
        if remove_repeating_header:
            print("Repeating header metadata removal enabled - only unique IDs will be stored")
    
    if paired_end:
        print(f"Paired-end mode enabled (output mode: {paired_end_mode})")
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (streaming mode)")

    # Premake lookup table
    BYTE_LOOKUP = np.array([f'{i:03d}'.encode('ascii') for i in range(256)])
    
    # For multiple flowcells tracking
    flowcell_metadata_list = []  # List including things like common_metadata, start_index, end_index
    current_flowcell_metadata = None
    flowcell_start_index = 0
    total_sequences = 0
    buffer = b''
    
    print("Collecting sequences...")
    
    with open(fastq_path, 'rb', buffering=32*1024*1024) as infile, \
         open(output_path, 'wb', buffering=32*1024*1024) as outfile:
        
        # Write SRA accession first if provided
        if sra_accession:
            outfile.write(f"@{sra_accession}\n".encode('utf-8'))
        
        # Reserve space for metadata headers (we'll come back to write these)
        metadata_position = outfile.tell()
        
        # Process file in chunks
        while True:
            chunk = infile.read(chunk_size_bytes)
            if not chunk and not buffer:
                break
            
            buffer += chunk
            
            # Find last complete record boundary
            last_at = buffer.rfind(b'\n@')
            if last_at == -1 or not chunk:  # Last chunk/no complete record
                process_buffer = buffer
                buffer = b''
            else:
                process_buffer = buffer[:last_at+1]
                buffer = buffer[last_at+1:]
            
            # Parse records from this buffer
            records, leftover, metadata, count = parse_fastq_records_from_buffer(
                process_buffer, total_sequences, base_map, phred_map,
                compress_headers, sequencer_type, paired_end,
                keep_bases, binary_bases, keep_quality, binary_quality
            )
            
            # Track flowcell changes if multiple_flowcells is enabled
            if metadata:
                if multiple_flowcells:
                    if current_flowcell_metadata is None:
                        # First flowcell
                        current_flowcell_metadata = metadata
                        flowcell_start_index = total_sequences
                    elif not metadata_dict_equals(current_flowcell_metadata, metadata):
                        # Flowcell changed! Save the previous flowcell range
                        flowcell_metadata_list.append((
                            current_flowcell_metadata.copy(),
                            flowcell_start_index,
                            total_sequences - 1
                        ))
                        current_flowcell_metadata = metadata
                        flowcell_start_index = total_sequences
                        print(f"Flowcell change detected at sequence {total_sequences}")
                elif multiple_flowcells and not metadata and current_flowcell_metadata is not None:
                    # Empty metadata but we had valid metadata before? So we'll skip this sequence/handle as error
                    print(f"WARNING: Failed to parse header at sequence {total_sequences}, skipping flowcell check")
                elif current_flowcell_metadata is None:
                    # Single flowcell mode (original behavior)
                    current_flowcell_metadata = metadata
            
            # Process and write records immediately
            if records:
                process_and_write_records(
                    records, outfile, base_map, quality_scaling,
                    custom_formula, phred_alphabet_max, min_quality, keep_bases,
                    binary_bases, binary, keep_quality,
                    remove_repeating_header, BYTE_LOOKUP
                )
            
            total_sequences += count
            buffer = leftover + buffer  # Keep leftover for next iteration
            
            if total_sequences % 100000 == 0:
                print(f"Collected {total_sequences:,} sequences...")
        
        # Finalize the last flowcell range if multiple flowcells enabled
        if multiple_flowcells and current_flowcell_metadata is not None:
            flowcell_metadata_list.append((
                current_flowcell_metadata.copy(),
                flowcell_start_index,
                total_sequences - 1
            ))
        elif current_flowcell_metadata:
            flowcell_metadata_list.append((current_flowcell_metadata, 0, total_sequences - 1))
        
        # Print flowcell summary
        if multiple_flowcells and len(flowcell_metadata_list) > 1:
            print(f"\nDetected {len(flowcell_metadata_list)} flowcells:")
            for idx, (metadata, start, end) in enumerate(flowcell_metadata_list):
                fc_id = metadata.get('flowcell', metadata.get('movie', 'unknown'))
                print(f"  Flowcell {idx + 1}: {fc_id} (sequences {start}-{end}, total: {end - start + 1})")
        
        # Now go back and write metadata headers at the beginning
        # Write metadata header(s) at the top if compression is enabled
        if compress_headers and sequencer_type != 'none':
            end_position = outfile.tell()
            outfile.seek(metadata_position)
            
            if multiple_flowcells and len(flowcell_metadata_list) > 0:
                # Write all flowcell metadata blocks
                for metadata, start_idx, end_idx in flowcell_metadata_list:
                    metadata_line = format_metadata_header(metadata, sequencer_type)
                    outfile.write(f"@{metadata_line}\n".encode('utf-8'))
                    
                    # Write parameter line with equation
                    equation = get_scaling_equation(quality_scaling, dataset_min_q, dataset_max_q, custom_formula, phred_alphabet_max)
                    parameter_line = f"#{equation}\n"
                    outfile.write(parameter_line.encode('utf-8'))
                    
                    # Write range line for multiple flowcells 
                    range_line = f"@RANGE:{start_idx}-{end_idx}\n"
                    outfile.write(range_line.encode('utf-8'))
            else:
                # Single flowcell (original behavior)
                if current_flowcell_metadata:
                    metadata_line = format_metadata_header(current_flowcell_metadata, sequencer_type)
                    outfile.write(f"@{metadata_line}\n".encode('utf-8'))
                
                # Write parameter line with equation
                equation = get_scaling_equation(quality_scaling, dataset_min_q, dataset_max_q, custom_formula, phred_alphabet_max)
                parameter_line = f"#{equation}\n"
                outfile.write(parameter_line.encode('utf-8'))
            
            # Fill remaining space to avoid corrupting data
            current_pos = outfile.tell()
            if current_pos < end_position:
                outfile.write(b'\n' * (end_position - current_pos))
    
    print(f"Total sequences written: {total_sequences:,}")

def main():
    argument_parser = argparse.ArgumentParser(description="Convert and compress FASTQ/FASTA files to scalar format")
    argument_parser.add_argument("input_path", type=str, help="Path of .fasta or .fastq file")
    argument_parser.add_argument("output_path", type=str, help="Output file path")

    # Quick mode implementation
    argument_parser.add_argument("mode", type=int,
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
                                choices=['none', 'illumina', 'pacbio', 'ont'],
                                default='none',
                                help="Sequencer type for header compression (default: none)")
    argument_parser.add_argument("--multiple_flowcells", type=int, default=0,
                                help="Enable multiple flowcell detection and tracking (0/1, default 0)")
    
    # Quality arguments
    argument_parser.add_argument("--extract_quality", type=int, default=1, 
                                help="For FASTQ: extract quality scores (0/1, default 1)")
    argument_parser.add_argument("--phred_offset", type=int, default=0,
                                help="Phred quality offset (if needed)(default 0)")
    argument_parser.add_argument("--min_quality", type=int, default=0,
                                help="Minimum quality score threshold (default 0)")
    argument_parser.add_argument("--quality_scaling", type=str, 
                                choices=['log','log_custom','log_adaptive', 'custom'], 
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
    argument_parser.add_argument("--phred_alphabet", type=str, default='phred42', # Added to prepare for future implementation of sequencing machine parameter
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
    export_scalars_to_txt(user_arguments.input_path, numpy_base_map, 
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
                          paired_end_mode=user_arguments.paired_end_mode)
    
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