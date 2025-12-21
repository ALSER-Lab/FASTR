import os
import numpy as np
import time
import argparse
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from multiprocessing import Pool
from functools import partial

@dataclass
class MetadataBlock:
    """Stores metadata for a flowcell/run at top of file"""
    common_metadata: Dict
    sequencer_type: str
    scaling_equation: str
    start_index: int
    end_index: int

def parse_metadata_header(data: bytes) -> Tuple[List[MetadataBlock], int, Optional[str]]:
    """
    Parse metadata headers from the beginning of the file in order to reconstruct it later
    """
    metadata_blocks = []
    sra_accession = None
    
    first_seq_marker = data.find(b'<')
    if first_seq_marker == -1:
        return metadata_blocks, 0, sra_accession
    
    # Find the @ before the first < (where sequence actually starts)
    search_start = max(0, first_seq_marker - 1000)  # Look back up to 1000 bytes
    header_before_seq = data[search_start:first_seq_marker]
    last_at = header_before_seq.rfind(b'\n@')
    
    if last_at != -1:
        # Found a header, sequence data starts there
        actual_data_start = search_start + last_at + 1  # +1 to skip the \n
    else:
        # No header found, might be at the very start
        if data[0:1] == b'@':
            actual_data_start = 0
        else:
            actual_data_start = first_seq_marker
    
    header_section = data[:actual_data_start].decode('utf-8', errors='ignore')
    lines = [line.strip() for line in header_section.split('\n') if line.strip()]
    
    print(f"Parsed {len(lines)} header lines (before first '<')")
    for i, line in enumerate(lines):
        print(f"  Line {i}: {line[:80]}")
    
    line_idx = 0
    
    # Check for SRA accession (first line starting with # but not a metadata tag)
    if lines and lines[0].startswith('#') and not any(x in lines[0] for x in ['SEQUENCER', 'RANGE', 'COMMON:']):
        # First line is SRA accession
        sra_accession = lines[0][1:]  # Remove #
        print(f"Found SRA accession: {sra_accession}")
        line_idx = 1
    
    # Only process lines that are actual metadata
    while line_idx < len(lines):
        line = lines[line_idx]
        
        # Stop if we hit sequence data (lines starting with @ but not @RANGE)
        if line.startswith('@') and not line.startswith('@RANGE'):
            # Check if this looks like a sequence header vs metadata
            test_content = line[1:]  # Remove @
            
            # If it's just numbers and colons/slashes (like "1:1" or "2/1"), it's a sequence header
            if re.match(r'^\d+[:/]\d+$', test_content):
                print(f"Detected sequence header at line {line_idx}: {line}")
                break
            
            # Otherwise we'll treat as metadata (shouldn't happen, just being safe though)
            print(f"WARNING: Found @ line that's not @RANGE: {line}")
            break
        
        # Stop if we somehow hit sequence data markers
        if line.startswith('<'):
            break
        
        # Look for common metadata line (new format)
        if line.startswith('#COMMON:'):
            common_metadata_str = line.split(':', 1)[1].strip()
            print(f"Found COMMON metadata: {common_metadata_str}")
            # We'll use this when we find the SEQUENCER line
            line_idx += 1
            continue
        
        # Look for sequencer type line
        if line.startswith('#SEQUENCER:'):
            sequencer_type = line.split(':', 1)[1].strip()
            
            # Look back one line for COMMON metadata
            common_metadata_str = None
            if line_idx > 0 and lines[line_idx - 1].startswith('#COMMON:'):
                common_metadata_str = lines[line_idx - 1].split(':', 1)[1].strip()
            
            line_idx += 1
            
            # Parse scaling equation (next line should start with #)
            scaling_equation = 'x'
            if line_idx < len(lines) and lines[line_idx].startswith('#') and not lines[line_idx].startswith('#SEQUENCER') and not lines[line_idx].startswith('#COMMON:'):
                scaling_equation = lines[line_idx][1:].strip()
                print(f"Found equation: {scaling_equation}")
                line_idx += 1
            
            # Parse common metadata based on sequencer type
            if common_metadata_str:
                common_metadata = parse_common_metadata(common_metadata_str, sequencer_type)
                print(f"Parsed common metadata: {common_metadata}")
            elif sra_accession:
                common_metadata = parse_common_metadata(sra_accession, 'srr')
            else:
                common_metadata = {}
            
            start_index = 0
            end_index = -1
            
            # Check for range
            if line_idx < len(lines) and lines[line_idx].startswith('@RANGE:'):
                range_str = lines[line_idx][7:]
                start_index, end_index = map(int, range_str.split('-'))
                print(f"Found range: {start_index}-{end_index}")
                line_idx += 1
            
            metadata_blocks.append(MetadataBlock(
                common_metadata=common_metadata,
                sequencer_type=sequencer_type,
                scaling_equation=scaling_equation,
                start_index=start_index,
                end_index=end_index
            ))
        else:
            line_idx += 1 # Skip any other lines
    
    return metadata_blocks, actual_data_start, sra_accession

def parse_custom_formula(formula: str, quality_scores: np.ndarray) -> np.ndarray:
    # Just copy-pasted from quality_processing.py 
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

def build_formula_func(formula: str):
    """Return a function f(x) that applies the custom formula."""
    cleaned = re.sub(r'^\s*f\s*\(\s*x\s*\)\s*=\s*', '', formula.strip()).replace('^', '**') # Normalize formula string once

    safe_dict = {
        'ln': np.log, 'log': np.log, 'log10': np.log10, 'exp': np.exp,
        'sqrt': np.sqrt, 'abs': np.abs, 'min': np.minimum, 'max': np.maximum,
        'np': np, '__builtins__': {}
    }

    def formula_func(x):
        local = dict(safe_dict)
        local['x'] = x
        with np.errstate(divide='ignore', invalid='ignore'):
            return eval(cleaned, local)

    return formula_func

def parse_common_metadata(metadata_str: str, sequencer_type: str) -> Dict:
    """Parse common metadata string based on sequencer type"""
    if sequencer_type == 'illumina':
        parts = metadata_str.split(':')
        return {
            'instrument': parts[0] if len(parts) > 0 else '',
            'run_id': parts[1] if len(parts) > 1 else '',
            'flowcell': parts[2] if len(parts) > 2 else '',
            'lane': parts[3] if len(parts) > 3 else ''
        }
    elif sequencer_type in ['pacbio_ccs', 'pacbio_hifi', 'pacbio_subread', 'pacbio_clr']:
        # Handle PacBio metadata: "movie:read_type" or just "movie"
        parts = metadata_str.split(':')
        result = {'movie': parts[0] if len(parts) > 0 else metadata_str}
        if len(parts) > 1:
            result['read_type'] = parts[1]
        return result
    elif sequencer_type == 'ont':
        kvs = {}
        for part in metadata_str.split(':'):
            if '=' in part:
                key, value = part.split('=', 1)
                kvs[key] = value
        return kvs
    elif sequencer_type == 'srr':
        match = re.match(r'([A-Z]+)(\d+)', metadata_str)
        if match:
            return {'prefix': match.group(1), 'accession': match.group(2)}
        return {}
    else:
        return {}


def reconstruct_header(unique_id: str, common_metadata: Dict, sequencer_type: str, 
                      pair_number: int = 0, sra_accession: Optional[str] = None) -> str:
    """Reconstruct full FASTQ header from FASTR"""
    if sequencer_type == 'illumina':
        # unique_id format: tile:x:y:read:filtered:control:barcode
        parts = unique_id.split(':')
        if len(parts) >= 7:
            header = f"@{common_metadata['instrument']}:{common_metadata['run_id']}:{common_metadata['flowcell']}:{common_metadata['lane']}:{parts[0]}:{parts[1]}:{parts[2]} {parts[3]}:{parts[4]}:{parts[5]}:{parts[6]}"
        else:
            header = f"@{unique_id}"
    
    elif sequencer_type == 'pacbio_ccs':
        # Format: @movie/zmw/ccs
        movie = common_metadata.get('movie', '')
        header = f"@{movie}/{unique_id}/ccs"
    
    elif sequencer_type == 'pacbio_hifi':
        # Format: @movie/zmw/ccs or @movie/zmw/ccs/fwd etc
        movie = common_metadata.get('movie', '')
        # unique_id might already include suffix like "123/fwd"
        if '/ccs' in unique_id:
            header = f"@{movie}/{unique_id}"
        else:
            header = f"@{movie}/{unique_id}/ccs"
    
    elif sequencer_type == 'pacbio_subread':
        # Format: @movie/zmw/start_end
        movie = common_metadata.get('movie', '')
        header = f"@{movie}/{unique_id}"
    
    elif sequencer_type == 'pacbio_clr':
        # Format: @movie/zmw/start_end
        movie = common_metadata.get('movie', '')
        header = f"@{movie}/{unique_id}"
    
    elif sequencer_type == 'ont':
        # unique_id format: read_id:key=value:key=value etc
        parts = unique_id.split(':', 1)
        read_id = parts[0]
        
        # Reconstruct w/ common metadata
        metadata_parts = []
        for key, value in common_metadata.items():
            metadata_parts.append(f"{key}={value}")
        
        if len(parts) > 1:
            metadata_parts.append(parts[1])
        
        header = f"@{read_id} {' '.join(metadata_parts)}"
    
    elif sequencer_type == 'srr':
        # unique_id format: read_index:spot
        parts = unique_id.split(':')
        if len(parts) >= 2:
            header = f"@{common_metadata['prefix']}{common_metadata['accession']}.{parts[0]} {parts[1]}"
        else:
            header = f"@{unique_id}"
    
    else:
        # If we get type of 'none', check if we have SRA accession
        if sra_accession:
            # Parse unique_id as read_num:pair_num
            parts = unique_id.split(':')
            if len(parts) >= 2:
                header = f"@{sra_accession}.{parts[0]} {parts[1]}"
            else:
                header = f"@{sra_accession}.{unique_id}"
        else:
            header = f"@{unique_id}"
    
    if pair_number > 0 and f"/{pair_number}" not in header and f" {pair_number}" not in header:
        header += f"/{pair_number}"
    
    return header


def reverse_base_map(gray_N=1, gray_A=63, gray_T=127, gray_C=191, gray_G=255) -> np.ndarray:
    """Create reverse mapping array from binary values to bases"""
    reverse_map = np.full(256, ord('N'), dtype=np.uint8)
    
    # map exact base vals
    reverse_map[gray_N] = ord('N')
    reverse_map[gray_A] = ord('A')
    reverse_map[gray_T] = ord('T')
    reverse_map[gray_C] = ord('C')
    reverse_map[gray_G] = ord('G')
    
    # When we have quality-scaled values, we need to map them to their initial base range
    # A: 1-63
    for i in range(1, 63):
        reverse_map[i] = ord('A')
    
    # T: 65-127
    for i in range(65, 127):
        reverse_map[i] = ord('T')

    # C: 129-191
    for i in range(129, 192):
        reverse_map[i] = ord('C')
    
    # G: 193-255
    for i in range(193, 255):
        reverse_map[i] = ord('G')
    
    return reverse_map

def reverse_scaling_to_quality(binary_values: np.ndarray,
                               base_ranges: dict, reverse_map: np.ndarray,
                               inverse_table: np.ndarray,
                               max_phred: int) -> np.ndarray:
    """
    Reverse the scaled quality values from the 63-per-base scheme.
    Returns reconstructed quality scores (0..max_phred).
    """
    
    # Map binary values to bases
    bases = reverse_map[binary_values]
    
    # Calculate y (scaled values 0-63) by subtracting base minimum
    y = np.zeros_like(binary_values, dtype=np.int32)
    for base_char, (lo, hi) in base_ranges.items():
        mask = bases == ord(base_char)
        if np.any(mask):
            y[mask] = binary_values[mask] - lo
    
    # Clamp y is within valid range
    # Weirdly enough np.clip clamps, not clips?
    y = np.clip(y, 0, 63)
    reconstructed_q = inverse_table[y] # Reconstruct quality scores w/ inverse lookup
    reconstructed_q = np.clip(reconstructed_q, 0, max_phred)
    
    return reconstructed_q

def build_inverse_quality_table(formula_func, max_phred, max_range=63):
    """Function designed to build the inverse table separately in order to not be looped within reverse_scaling_to_quality"""
    q_possible = np.arange(max_phred + 1, dtype=np.float32)
    
    # Apply forward formula to get scaled values
    scaled_for_q = formula_func(q_possible).astype(np.float32)
    scaled_for_q = np.nan_to_num(scaled_for_q, nan=0.0, posinf=max_range, neginf=0.0)
    scaled_for_q = np.clip(scaled_for_q, 0, max_range) # Clamp
    scaled_int_for_q = np.round(scaled_for_q).astype(int)  # Round to uint8
    inverse_table = np.full(max_range + 1, -1, dtype=int)
    
    # For each quality score, store it at its corresponding scaled position
    for q, s in zip(q_possible, scaled_int_for_q):
        if 0 <= s <= max_range:
            inverse_table[s] = int(q)
    
    valid_indices = np.where(inverse_table != -1)[0] # Fill gaps in inverse table 
    
    if len(valid_indices) > 0:
        # Fill beginning
        if valid_indices[0] > 0:
            inverse_table[:valid_indices[0]] = inverse_table[valid_indices[0]]
        
        # Fill gaps between valid values
        for i in range(len(valid_indices) - 1):
            start_idx = valid_indices[i]
            end_idx = valid_indices[i + 1]
            if end_idx - start_idx > 1:
                start_val = inverse_table[start_idx]
                end_val = inverse_table[end_idx]
                gap_size = end_idx - start_idx
                for j in range(1, gap_size):
                    interpolated = start_val + (end_val - start_val) * j / gap_size
                    inverse_table[start_idx + j] = int(np.round(interpolated))
        
        # Fill end
        if valid_indices[-1] < max_range:
            inverse_table[valid_indices[-1] + 1:] = inverse_table[valid_indices[-1]]
    else:
        # If no valid mappings found, use default quality of 0
        inverse_table[:] = 0
    
    return inverse_table

def quality_to_ascii(quality_scores: np.ndarray, phred_offset: int = 33) -> bytes:
    """Convert numeric quality scores to ASCII string"""
    return bytes((quality_scores + phred_offset).astype(np.uint8))


def process_chunk_worker_reconstruction(chunk_data, reverse_map, base_ranges, 
                                       metadata_blocks, inverse_tables, 
                                       phred_alphabet_max, phred_offset, sra_accession):
    """
    Worker function that processes a single chunk of binary sequence data in parallel.
    We use same streaming + multiprocessing approach as in our to_fastr converison methodology (can be more specifically found in chunk_processor.py)
    (chunk_id, reconstructed_fastq_text, sequence_count)
    """
    try:
        chunk_id, chunk_binary, start_seq_idx = chunk_data
        print(f"Worker processing chunk {chunk_id} starting at sequence {start_seq_idx}")
        
        output_buffer = []
        sequence_count = start_seq_idx
        pos = 0
        
        # Determine which metadata block to use for this chunk
        # (for files with multiple flowcells, different ranges use different metadata)
        current_metadata_idx = 0
        for i, mb in enumerate(metadata_blocks):
            if mb.end_index == -1 or sequence_count <= mb.end_index:
                current_metadata_idx = i
                break
        
        current_metadata = metadata_blocks[current_metadata_idx]
        current_inverse_table = inverse_tables[current_metadata_idx]
        
        # Process all sequences in this chunk
        while pos < len(chunk_binary):
            # Find sequence start marker '<'
            if chunk_binary[pos:pos+1] != b'<':
                pos += 1
                continue
            
            pos += 1  # Skip <
            
            # Find newline to determine sequence length
            seq_end = chunk_binary.find(b'\n', pos)
            if seq_end == -1:
                if pos < len(chunk_binary):
                    seq_end = len(chunk_binary)
                else:
                    break
            
            # Extract binary sequence data
            seq_data = chunk_binary[pos:seq_end]
            
            # Check if we need to switch metadata blocks (for multiple flowcells)
            if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                next_metadata = metadata_blocks[current_metadata_idx + 1]
                if sequence_count >= next_metadata.start_index:
                    current_metadata = next_metadata
                    current_inverse_table = inverse_tables[current_metadata_idx + 1]
                    current_metadata_idx += 1
            
            # Convert binary values back to bases
            seq_array = np.frombuffer(seq_data, dtype=np.uint8)
            bases_array = reverse_map[seq_array]
            bases = bases_array.tobytes().decode('ascii')
            
            # Try to find header before this sequence by searching backwards
            header_search_start = max(0, pos - 500)
            header_section = chunk_binary[header_search_start:pos]

            unique_id = None
            pair_number = 0

            # Find last '@' before '<' (the header for this sequence)
            last_at = header_section.rfind(b'@')
            if last_at != -1:
                header_line_start = header_search_start + last_at
                header_line_end = chunk_binary.find(b'\n', header_line_start)
                if header_line_end != -1 and header_line_end <= pos:
                    header_content = chunk_binary[header_line_start+1:header_line_end].decode('utf-8', errors='ignore').strip()
                    
                    # Check for pair number suffix (for paired-end reads)
                    if '/' in header_content:
                        parts = header_content.rsplit('/', 1)
                        unique_id = parts[0]
                        try:
                            pair_number = int(parts[1])
                        except:
                            unique_id = header_content
                    else:
                        unique_id = header_content

            # A bug I ran into when implementing multiprocessing was the first sequence header being improperly reconstructed (normally as "seq_0")
            # A duct-tape-esque fix right now is that if we're very near the start and haven't found a header yet, search from beginning
            elif pos <= 510 and chunk_binary[:pos].count(b'<') == 1:
                first_at = chunk_binary.find(b'@')
                if first_at != -1 and first_at < pos:
                    header_line_end = chunk_binary.find(b'\n', first_at)
                    print(f"DEBUG: header_line_end = {header_line_end}")
                    if header_line_end != -1 and header_line_end < pos:
                        header_content = chunk_binary[first_at+1:header_line_end].decode('utf-8', errors='ignore').strip()
                        
                        # Check for pair number
                        if '/' in header_content:
                            parts = header_content.rsplit('/', 1)
                            unique_id = parts[0]
                            try:
                                pair_number = int(parts[1])
                            except:
                                unique_id = header_content
                        else:
                            unique_id = header_content
            
            # Reconstruct full header from compressed unique id
            if unique_id:
                header = reconstruct_header(unique_id, current_metadata.common_metadata,
                                          current_metadata.sequencer_type, pair_number, sra_accession)
            else:
                # Fallback is togenerate generic header if we couldn't find the compressed one
                if sra_accession:
                    header = f"@{sra_accession}.{sequence_count}"
                else:
                    header = f"@seq_{sequence_count}"
            
            # Apply reverse scaling to reconstruct original quality scores
            quality_scores = reverse_scaling_to_quality(
                seq_array,
                base_ranges,
                reverse_map,
                current_inverse_table,
                max_phred=phred_alphabet_max
            )
            quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
            
            # Add complete FASTQ record to output buffer
            output_buffer.append(f"{header}\n{bases}\n+\n{quality_string}\n")
            
            sequence_count += 1
            pos = seq_end + 1
        
        print(f"Worker completed chunk {chunk_id}, processed {sequence_count - start_seq_idx} sequences")
        return (chunk_id, ''.join(output_buffer), sequence_count - start_seq_idx)
    
    except Exception as e:
        print(f"ERROR in worker processing chunk {chunk_id}: {e}")
        raise

def reconstruct_fastq(input_path: str, output_path: str, 
                     gray_N: int = 1, gray_A: int = 63, gray_T: int = 127,
                     gray_C: int = 191, gray_G: int = 255, phred_alphabet_max: int = 41,
                     phred_offset: int = 33, chunk_size_mb: int = 32, num_workers: int = 4):
    """
    Reconstruct FASTQ file from FASTR (mode 3) using parallel processing.
    Assumes binary writing was used in sequence_converter (--binary_write 1).
    """
    print(f"Reading FASTR: {input_path}")
    
    # Read entire file as binary
    with open(input_path, 'rb') as f:
        data = f.read()
    file_size = len(data)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**2):.2f} MB)")
    print(f"Using {num_workers} worker processes for parallel reconstruction")
    
    # Parse metadata headers from the beginning of the file
    metadata_blocks, data_start_byte, sra_accession = parse_metadata_header(data)
    
    if not metadata_blocks:
        print("ERROR: No metadata found?")
        return
    
    print(f"Found {len(metadata_blocks)} metadata block(s)")
    if sra_accession:
        print(f"SRA Accession: {sra_accession}")
    
    for i, mb in enumerate(metadata_blocks):
        print(f"  Block {i+1}: {mb.sequencer_type}, equation: {mb.scaling_equation}")
    
    # Create reverse base mapping (binary values:ASCII bases)
    reverse_map = reverse_base_map(gray_N, gray_A, gray_T, gray_C, gray_G)
    
    # Pre-compute inverse quality tables for each metadata block (making it easier to reverse)
    inverse_tables = []
    for mb in metadata_blocks:
        formula_func = build_formula_func(mb.scaling_equation)
        inverse_table = build_inverse_quality_table(formula_func, phred_alphabet_max)
        inverse_tables.append(inverse_table)
    
    # Base ranges for quality decoding (each base has a 63-value range)
    base_ranges = {
        'A': (1, 63),
        'T': (65, 127),
        'C': (129, 191),
        'G': (193, 255),
        'N': (1, 1)
    }
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    # Generator function to yield chunks from binary data
    def chunk_generator():
        buffer = data[data_start_byte:]
        chunk_id = 0
        start_seq_idx = 0
        pos = 0
        
        while pos < len(buffer):
            chunk_end = min(pos + chunk_size_bytes, len(buffer))
            
            # Find last complete sequence boundary (look for < followed by content and newline)
            # We want to split on sequence boundaries to avoid cutting sequences in half
            search_start = max(pos, chunk_end - 1000)
            last_marker = buffer.rfind(b'<', search_start, chunk_end)
            
            if last_marker != -1 and last_marker > pos:
                # Find the newline after this sequence to get the complete sequence
                seq_end = buffer.find(b'\n', last_marker)
                if seq_end != -1 and seq_end < len(buffer):
                    chunk_end = seq_end + 1
            
            # Include overlap from before the chunk to capture headers
            # We do this because headers are stored before their sequences, 
            # so we need to look back to check whether each chunk contains the headers for all its sequences
            overlap_start = max(0, pos - 500)
            chunk_binary = buffer[overlap_start:chunk_end]
            
            if chunk_binary:
                # Estimate sequences in this chunk for tracking purposes
                # Only count sequences in the NEW data (not the overlap) to avoid double counting
                actual_start_in_chunk = pos - overlap_start
                estimated_seqs = chunk_binary[actual_start_in_chunk:].count(b'<')
                
                yield (chunk_id, chunk_binary, start_seq_idx)
                
                start_seq_idx += estimated_seqs
                chunk_id += 1
                
                if chunk_id % 10 == 0:
                    print(f"Read {chunk_id} chunks...")
            
            pos = chunk_end
            
            if pos >= len(buffer):
                break
    
    print("Reconstructing FASTQ with parallel processing...")
    total_sequences = 0
    
    with open(output_path, 'w', buffering=chunk_size_bytes) as outfile:
        # Create worker pool and process chunks as they come in
        with Pool(processes=num_workers) as pool:
            # Create partial function with fixed arguments
            worker_func = partial(
                process_chunk_worker_reconstruction,
                reverse_map=reverse_map,
                base_ranges=base_ranges,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                sra_accession=sra_accession
            )
            
            # Use imap to process chunks as they're read (streaming), like how it is done in to_fastr.py
            # chunksize=1 makes order preserved for sequential writing
            for chunk_id, fastq_text, count in pool.imap(worker_func, chunk_generator(), chunksize=1):
                # Write immediately as each chunk completes
                outfile.write(fastq_text)
                
                total_sequences += count
                
                if total_sequences % 100000 == 0:
                    print(f"Reconstructed {total_sequences:,} sequences...")
    
    print(f"\nTotal sequences reconstructed: {total_sequences:,}")
    print(f"Output saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Reconstruct FASTQ files from compressed binary format (mode 3)")
    parser.add_argument("input_path", type=str, help="Path to compressed binary file")
    parser.add_argument("output_path", type=str, help="Output FASTQ file path")
    
    # Quality parameters
    parser.add_argument("--phred_offset", type=int, default=33,
                        help="Phred quality offset for output (default: 33)")
    parser.add_argument("--phred_alphabet", type=str, default='phred42',
                       help="Phred alphabet used (phred42/phred63/phred94)")
    
    # Base mapping (copy-pasted from sequence_converter.py)
    parser.add_argument("--gray_N", type=int, default=1)
    parser.add_argument("--gray_A", type=int, default=63)
    parser.add_argument("--gray_T", type=int, default=127)
    parser.add_argument("--gray_C", type=int, default=191)
    parser.add_argument("--gray_G", type=int, default=255)
    
    # Multiprocessing arguments
    parser.add_argument("--chunk_size_mb", type=int, default=32,
                        help="Chunk size in MB for parallel processing (default: 32)")
    parser.add_argument("--num_workers", type=int, default=4,
                        help="Number of parallel workers (default: 4, use 4+ for large files)")
    
    args = parser.parse_args()
    
    # Determine phred alphabet max
    # This may not properly be transferred to our reverse_scaling_quality(), but it is fine for now. 
    if args.phred_alphabet == "phred42":
        phred_alphabet_max = 41
    if args.phred_alphabet == "phred63":
        phred_alphabet_max = 62
    elif args.phred_alphabet == "phred94":
        phred_alphabet_max = 93
    
    start_time = time.perf_counter()
    
    reconstruct_fastq(
        args.input_path, args.output_path,
        gray_N=args.gray_N, gray_A=args.gray_A,
        gray_T=args.gray_T, gray_C=args.gray_C,
        gray_G=args.gray_G,
        phred_alphabet_max=phred_alphabet_max,
        phred_offset=args.phred_offset,
        chunk_size_mb=args.chunk_size_mb,
        num_workers=args.num_workers
    )
    
    end_time = time.perf_counter()
    print(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()