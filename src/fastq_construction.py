import os
import numpy as np
import time
import argparse
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from multiprocessing import Pool
from functools import partial
import mmap

@dataclass
class MetadataBlock:
    """Stores metadata for a flowcell/run at top of file"""
    structure_template: str
    sequencer_type: str
    scaling_equation: str
    start_index: int
    end_index: int

def find_structure_prefix(structure_template: str) -> str:
    """
    Extract the constant prefix from structure template before first {REPEATING_X}.
    """
    if not structure_template:
        return ""
    
    # Find first occurrence of {REPEATING_
    match = re.search(r'\{REPEATING_\d+\}', structure_template)
    
    if match:
        prefix = structure_template[:match.start()] # Return everything before the first placeholder
        if prefix.endswith(':') or prefix.endswith('/'): # Remove trailing delimiter if present
            prefix = prefix[:-1]
        return prefix
    
    return structure_template


def parse_metadata_header(data: bytes, mode: int) -> Tuple[List[MetadataBlock], int, Optional[str]]:
    """
    Parse metadata headers from the beginning of the file in order to reconstruct it later
    """
    metadata_blocks = []
    sra_accession = None
    
    if mode == 1:
        # Mode 1: Only bases compressed (no header compression)
        # Parse metadata for quality reconstruction, then find where actual sequences start

        # Fixed bug that was caused by mode 1 "not needing to look at the metadata" despite it needing to reverse the scaling equation
        # All modes look at the metadata header
        first_header = data.find(b'\n@')
        if first_header == -1:
            first_header = data.find(b'@')
            if first_header == -1:
                first_header = 0
        
        # Parse metadata section before first @ header
        if first_header > 0:
            header_section = data[:first_header].decode('utf-8', errors='ignore')
            lines = [line.strip() for line in header_section.split('\n') if line.strip()]
            
            print(f"Parsed {len(lines)} header lines")
            for i, line in enumerate(lines[:5]):  # Show first 5 lines
                print(f"  Line {i}: {line[:80]}")
            
            line_idx = 0
            
            # Check for SRA accession
            if lines and lines[0].startswith('#') and not any(x in lines[0] for x in ['SEQUENCER', 'RANGE', 'STRUCTURE:']):
                sra_accession = lines[0][1:]
                print(f"Found SRA accession: {sra_accession}")
                line_idx = 1
            
            # Process metadata lines
            while line_idx < len(lines):
                line = lines[line_idx]
                
                if line.startswith('@'):
                    print(f"Reached sequence headers at line {line_idx}")
                    break
                
                if line.startswith('<'):
                    break
                
                if line.startswith('#STRUCTURE:'):
                    structure_template = line.split(':', 1)[1].strip()
                    print(f"Found STRUCTURE metadata: {structure_template}")
                    line_idx += 1
                    continue
                
                if line.startswith('#SEQUENCER:'):
                    sequencer_type = line.split(':', 1)[1].strip()
                    
                    structure_template = None
                    if line_idx > 0 and lines[line_idx - 1].startswith('#STRUCTURE:'):
                        structure_template = lines[line_idx - 1].split(':', 1)[1].strip()
                    
                    line_idx += 1
                    
                    scaling_equation = 'x'
                    if line_idx < len(lines) and lines[line_idx].startswith('#QUAL_SCALE:'):
                        scaling_equation = lines[line_idx].split(':', 1)[1].strip()
                        print(f"Found equation: {scaling_equation}")
                        line_idx += 1
                    
                    start_index = 0
                    end_index = -1
                    
                    if line_idx < len(lines) and lines[line_idx].startswith('#RANGE:'):
                        range_str = lines[line_idx].split(':', 1)[1].strip()
                        start_index, end_index = map(int, range_str.split('-'))
                        print(f"Found range: {start_index}-{end_index}")
                        line_idx += 1
                    
                    metadata_blocks.append(MetadataBlock(
                        structure_template=structure_template or "",
                        sequencer_type=sequencer_type,
                        scaling_equation=scaling_equation,
                        start_index=start_index,
                        end_index=end_index
                    ))
                else:
                    line_idx += 1
        
        # Find actual data start (first @ header)
        actual_data_start = first_header if first_header > 0 else 0
        return metadata_blocks, actual_data_start, sra_accession
    
    # Mode 0, 2, and 3: Header compression enabled
    first_seq_marker = data.find(b'<') if mode in [2, 3] else data.find(b'\n@')
    if first_seq_marker == -1:
        return metadata_blocks, 0, sra_accession
    
    # Find the @ before the first < (where sequence actually starts)
    search_start = max(0, first_seq_marker - 1000)
    header_before_seq = data[search_start:first_seq_marker]
    last_at = header_before_seq.rfind(b'\n@')
    
    if last_at != -1:
        actual_data_start = search_start + last_at + 1
    else:
        if data[0:1] == b'@':
            actual_data_start = 0
        else:
            actual_data_start = first_seq_marker
    
    header_section = data[:actual_data_start].decode('utf-8', errors='ignore')
    lines = [line.strip() for line in header_section.split('\n') if line.strip()]
    
    print(f"Parsed {len(lines)} header lines")
    for i, line in enumerate(lines[:5]):  # Show first 5 lines
        print(f"  Line {i}: {line[:80]}")
    
    line_idx = 0
    
    # Check for SRA accession
    if lines and lines[0].startswith('#') and not any(x in lines[0] for x in ['SEQUENCER', 'RANGE', 'STRUCTURE:']):
        sra_accession = lines[0][1:]
        print(f"Found SRA accession: {sra_accession}")
        line_idx = 1
    
    # Process metadata lines
    while line_idx < len(lines):
        line = lines[line_idx]
        
        if line.startswith('@') and not line.startswith('#RANGE:'):
            test_content = line[1:]
            if re.match(r'^\d+[:/]\d+$', test_content):
                print(f"Detected sequence header at line {line_idx}: {line}")
                break
            print(f"WARNING: Found @ line that's not #RANGE: {line}")
            break
        
        if line.startswith('<'):
            break
        
        if line.startswith('#STRUCTURE:'):
            structure_template = line.split(':', 1)[1].strip()
            print(f"Found STRUCTURE metadata: {structure_template}")
            line_idx += 1
            continue
        
        if line.startswith('#SEQUENCER:'):
            sequencer_type = line.split(':', 1)[1].strip()
                    
            structure_template = None
            if line_idx > 0 and lines[line_idx - 1].startswith('#STRUCTURE:'):
                structure_template = lines[line_idx - 1].split(':', 1)[1].strip()
                line_idx += 1
                scaling_equation = 'x'
                if line_idx < len(lines) and lines[line_idx].startswith('#QUAL_SCALE:'):
                    scaling_equation = lines[line_idx].split(':', 1)[1].strip()
                    print(f"Found equation: {scaling_equation}")
                    line_idx += 1
            
            start_index = 0
            end_index = -1
            
            if line_idx < len(lines) and lines[line_idx].startswith('#RANGE:'):
                range_str = lines[line_idx].split(':', 1)[1].strip()
                start_index, end_index = map(int, range_str.split('-'))
                print(f"Found range: {start_index}-{end_index}")
                line_idx += 1
            
            metadata_blocks.append(MetadataBlock(
                structure_template=structure_template or "",
                sequencer_type=sequencer_type,
                scaling_equation=scaling_equation,
                start_index=start_index,
                end_index=end_index
            ))
        else:
            line_idx += 1
    
    return metadata_blocks, actual_data_start, sra_accession


def get_delimiter_for_sequencer(sequencer_type: str) -> str:
    """Get the delimiter character used by each sequencer type"""
    if sequencer_type == 'illumina' or sequencer_type == 'old_illumina':
        return ':'
    elif sequencer_type.startswith('pacbio'):
        return '/'
    elif sequencer_type == 'ont':
        return ':'
    elif sequencer_type == 'srr':
        return ':'
    return ':'


def reconstruct_header_from_structure(structure: str, unique_id: str, sequencer_type: str, pair_number: int = 0) -> str:
    """
    Reconstruct full header from structure template and unique ID.
    For adaptive format, detects delimiter from structure.
    """
    if sequencer_type == 'adaptive':
        # Detect primary delimiter from structure
        if ' ' in structure:
            primary_delimiter = ' '
        elif '/' in structure:
            primary_delimiter = '/'
        elif ':' in structure:
            primary_delimiter = ':'
        else:
            primary_delimiter = ' '
        
        unique_parts = unique_id.split(primary_delimiter)
        
        if primary_delimiter == ' ' and '.' in structure.split(' ')[0]:
            structure_temp = structure.replace('.', ' ')
            structure_parts = structure_temp.split(' ')
            
            result_parts = []
            unique_idx = 0
            for s_part in structure_parts:
                if s_part.startswith('{REPEATING_'):
                    if unique_idx < len(unique_parts):
                        result_parts.append(unique_parts[unique_idx])
                        unique_idx += 1
                else:
                    result_parts.append(s_part)
            
            if len(result_parts) >= 2:
                result = result_parts[0] + '.' + ' '.join(result_parts[1:])
            else:
                result = ' '.join(result_parts)
        else:
            # Normal single-delimiter reconstruction
            result = structure
            for i, part in enumerate(unique_parts, 1):
                placeholder = f"{{REPEATING_{i}}}"
                result = result.replace(placeholder, part, 1)
    else: # Normal sequencer types
        delimiter = get_delimiter_for_sequencer(sequencer_type)
        unique_parts = unique_id.split(delimiter)
        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part)
    
    if pair_number > 0:
        result = f"{result}/{pair_number}"
    
    return f"@{result}"


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


def create_base_map():
    base_table = np.zeros(256, dtype=np.int32)
    base_table[1:64] = 1
    base_table[65:128] = 65
    base_table[129:192] = 129
    base_table[193:256] = 193
    return base_table

def reverse_base_map(gray_N=1, gray_A=63, gray_T=127, gray_C=191, gray_G=255):
    reverse_map = np.full(256, ord('N'), dtype=np.uint8)
    reverse_map[1:64] = ord('A')
    reverse_map[65:128] = ord('T')
    reverse_map[129:193] = ord('C')
    reverse_map[193:256] = ord('G')
    reverse_map[[gray_N, gray_A, gray_T, gray_C, gray_G]] = [ord('N'), ord('A'), ord('T'), ord('C'), ord('G')]
    
    return reverse_map


def reverse_scaling_to_quality(binary_values: np.ndarray,
                                       subtract_table: np.ndarray,
                                       inverse_table: np.ndarray,
                                       max_phred: int) -> np.ndarray:
    """
    Reverse the scaled quality values from the 63-per-base scheme.
    Returns reconstructed quality scores (0..max_phred).
    """
    y = binary_values - subtract_table[binary_values]
    y = np.clip(y, 0, 63)
    
    # Single inverse lookup
    reconstructed_q = inverse_table[y]
    return np.clip(reconstructed_q, 0, max_phred)

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


def process_chunk_worker_reconstruction(chunk_data, reverse_map, subtract_table, 
                                       base_ranges, metadata_blocks, inverse_tables, 
                                       phred_alphabet_max, phred_offset, sra_accession,
                                       mode):
    """
    Worker function that processes a single chunk of binary sequence data in parallel.
    """
    try:
        chunk_id, chunk_binary, start_seq_idx = chunk_data
        print(f"Worker processing chunk {chunk_id} starting at sequence {start_seq_idx}")
        
        output_buffer = []
        sequence_count = start_seq_idx
        pos = 0
        
        # Determine which metadata block to use for this chunk
        current_metadata_idx = 0
        if metadata_blocks:
            for i, mb in enumerate(metadata_blocks):
                if mb.end_index == -1 or sequence_count <= mb.end_index:
                    current_metadata_idx = i
                    break
            current_metadata = metadata_blocks[current_metadata_idx]
            current_inverse_table = inverse_tables[current_metadata_idx] if inverse_tables else None
        else:
            current_metadata = None
            current_inverse_table = None
        
        # Mode 1: Only bases compressed (original headers intact)
        if mode == 1:
            lines = chunk_binary.split(b'\n')
            i = 0
            while i < len(lines):
                if lines[i].startswith(b'@'):
                    # Original header
                    header = lines[i].decode('utf-8', errors='ignore')
                    i += 1
                    if i < len(lines) and lines[i].startswith(b'<'):
                        # Binary sequence
                        seq_data = lines[i][1:]  # Skip '<'
                        seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                        bases_array = reverse_map[seq_array]
                        bases = bases_array.tobytes().decode('ascii')
                        
                        # Reconstruct quality if inverse table available
                        if current_inverse_table is not None:
                            quality_scores = reverse_scaling_to_quality(
                                seq_array, subtract_table,
                                current_inverse_table, max_phred=phred_alphabet_max
                            )
                            quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
                        else:
                            # Generate default quality scores
                            quality_string = 'I' * len(bases)
                        
                        output_buffer.append(f"{header}\n{bases}\n+\n{quality_string}\n")
                        sequence_count += 1
                        i += 1
                    else:
                        i += 1
                else:
                    i += 1
        
        # Mode 0: Only headers compressed (bases as text)
        elif mode == 0:
            lines = chunk_binary.split(b'\n')
            i = 0
            while i < len(lines):
                # Skip empty lines
                if not lines[i].strip():
                    i += 1
                    continue
                    
                if lines[i].startswith(b'@'):
                    # Compressed header - reconstruct it
                    header_content = lines[i][1:].decode('utf-8', errors='ignore').strip()
                    
                    pair_number = 0
                    if '/' in header_content:
                        parts = header_content.rsplit('/', 1)
                        unique_id = parts[0]
                        try:
                            pair_number = int(parts[1])
                        except:
                            unique_id = header_content
                    else:
                        unique_id = header_content
                    
                    if current_metadata and current_metadata.structure_template:
                        header = reconstruct_header_from_structure(
                            current_metadata.structure_template,
                            unique_id,
                            current_metadata.sequencer_type,
                            pair_number
                        )
                    else:
                        header = f"@{unique_id}"
                    
                    i += 1
                    # Get sequence line (skip if empty)
                    while i < len(lines) and not lines[i].strip():
                        i += 1
                    
                    if i < len(lines):
                        # Text bases
                        bases = lines[i].decode('ascii', errors='ignore').strip()
                        i += 1
                        
                        # Skip '+' line (and any empty lines before it)
                        while i < len(lines) and not lines[i].strip():
                            i += 1
                        if i < len(lines) and lines[i].strip() == b'+':
                            i += 1
                        
                        # Skip empty lines before quality
                        while i < len(lines) and not lines[i].strip():
                            i += 1
                            
                        # Get quality line
                        if i < len(lines):
                            quality_string = lines[i].decode('ascii', errors='ignore').strip()
                            output_buffer.append(f"{header}\n{bases}\n+\n{quality_string}\n".encode('ascii'))
                            sequence_count += 1
                            i += 1
                else:
                    i += 1

        # Mode 3: No headers at all, just <bases on each line
        elif mode == 3:
            lines = chunk_binary.split(b'\n')
            i = 0
            while i < len(lines):
                if lines[i].startswith(b'<'):
                    # Binary sequence
                    seq_data = lines[i][1:]  # Skip '<'
                    
                    # Check if we need to switch metadata blocks
                    if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                        next_metadata = metadata_blocks[current_metadata_idx + 1]
                        if sequence_count >= next_metadata.start_index:
                            current_metadata = next_metadata
                            current_inverse_table = inverse_tables[current_metadata_idx + 1]
                            current_metadata_idx += 1
                    
                    seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                    bases_array = reverse_map[seq_array]
                    bases = bases_array.tobytes().decode('ascii')
                    
                    # Generate header from structure template prefix only
                    if current_metadata and current_metadata.structure_template:
                        # Extract only the constant prefix before first {REPEATING_X}
                        header_prefix = find_structure_prefix(current_metadata.structure_template)
                        header = f"@{header_prefix}" if header_prefix else "@seq"
                    elif sra_accession:
                        header = f"@{sra_accession}"
                    else:
                        header = "@seq"
                    
                    # Reconstruct quality
                    if current_inverse_table is not None:
                        quality_scores = reverse_scaling_to_quality(
                            seq_array, subtract_table,
                            current_inverse_table, max_phred=phred_alphabet_max
                        )
                        quality_string = bytes((quality_scores + phred_offset).astype(np.uint8)).decode('ascii')
                    else:
                        quality_string = 'I' * len(bases)
                    
                    output_buffer.append(f"{header}\n".encode('ascii'))
                    output_buffer.append(bases.encode('ascii'))
                    output_buffer.append(b"\n+\n")
                    output_buffer.append(quality_string.encode('ascii'))
                    output_buffer.append(b"\n")
                    sequence_count += 1
                
                i += 1

        
        # Mode 2: Headers compressed, bases as binary (original behavior)
        else:
            while pos < len(chunk_binary):
                if chunk_binary[pos:pos+1] != b'<':
                    pos += 1
                    continue
                
                pos += 1
                
                seq_end = chunk_binary.find(b'\n', pos)
                if seq_end == -1:
                    if pos < len(chunk_binary):
                        seq_end = len(chunk_binary)
                    else:
                        break
                
                seq_data = chunk_binary[pos:seq_end]
                
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_inverse_table = inverse_tables[current_metadata_idx + 1]
                        current_metadata_idx += 1
                
                seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                bases_array = reverse_map[seq_array]
                bases = bases_array.tobytes().decode('ascii')
                
                header_search_start = max(0, pos - 500)
                header_section = chunk_binary[header_search_start:pos]

                unique_id = None
                pair_number = 0

                last_at = header_section.rfind(b'@')
                if last_at != -1:
                    header_line_start = header_search_start + last_at
                    header_line_end = chunk_binary.find(b'\n', header_line_start)
                    if header_line_end != -1 and header_line_end <= pos:
                        header_content = chunk_binary[header_line_start+1:header_line_end].decode('utf-8', errors='ignore').strip()
                        
                        if '/' in header_content:
                            parts = header_content.rsplit('/', 1)
                            unique_id = parts[0]
                            try:
                                pair_number = int(parts[1])
                            except:
                                unique_id = header_content
                        else:
                            unique_id = header_content

                elif pos <= 510 and chunk_binary[:pos].count(b'<') == 1:
                    first_at = chunk_binary.find(b'@')
                    if first_at != -1 and first_at < pos:
                        header_line_end = chunk_binary.find(b'\n', first_at)
                        if header_line_end != -1 and header_line_end < pos:
                            header_content = chunk_binary[first_at+1:header_line_end].decode('utf-8', errors='ignore').strip()
                            
                            if '/' in header_content:
                                parts = header_content.rsplit('/', 1)
                                unique_id = parts[0]
                                try:
                                    pair_number = int(parts[1])
                                except:
                                    unique_id = header_content
                            else:
                                unique_id = header_content
                
                if unique_id and current_metadata and current_metadata.structure_template:
                    header = reconstruct_header_from_structure(
                        current_metadata.structure_template,
                        unique_id,
                        current_metadata.sequencer_type,
                        pair_number
                    )
                else:
                    if sra_accession:
                        header = f"@{sra_accession}"
                    else:
                        header = f"@seq"
                
                if current_inverse_table is not None:
                    quality_scores = reverse_scaling_to_quality(
                        seq_array, subtract_table,
                        current_inverse_table, max_phred=phred_alphabet_max
                    )
                    quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
                else:
                    quality_string = 'I' * len(bases)
                
                output_buffer.append(header.encode('ascii'))
                output_buffer.append(b'\n')
                output_buffer.append(bases.encode('ascii'))
                output_buffer.append(b'\n+\n')
                output_buffer.append(quality_string.encode('ascii'))
                output_buffer.append(b'\n')
                                
                sequence_count += 1
                pos = seq_end + 1
            
        
        print(f"Worker completed chunk {chunk_id}, processed {sequence_count - start_seq_idx} sequences")
        return (chunk_id, b''.join(output_buffer), sequence_count - start_seq_idx)
    
    except Exception as e:
        print(f"ERROR in worker processing chunk {chunk_id}: {e}")
        import traceback
        traceback.print_exc()
        raise


def reconstruct_fastq(input_path: str, output_path: str, 
                     gray_N: int = 1, gray_A: int = 63, gray_T: int = 127,
                     gray_C: int = 191, gray_G: int = 255, phred_alphabet_max: int = 41,
                     phred_offset: int = 33, chunk_size_mb: int = 32, num_workers: int = 4,
                     mode: int = 2):
    """
    Reconstruct FASTQ file from FASTR using parallel processing.
    """
    print(f"Reading FASTR: {input_path}")
    print(f"Reconstruction mode: {mode}")
    
    with open(input_path, 'rb') as f:
        data = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    file_size = len(data)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**2):.2f} MB)")
    print(f"Using {num_workers} worker processes for parallel reconstruction")
    
    metadata_blocks, data_start_byte, sra_accession = parse_metadata_header(data, mode)
    
    if mode in [0, 2] and not metadata_blocks:
        print("WARNING: No metadata found for header reconstruction")
    
    if metadata_blocks:
        print(f"Found {len(metadata_blocks)} metadata block(s)")
        for i, mb in enumerate(metadata_blocks):
            print(f"  Block {i+1}: {mb.sequencer_type}, equation: {mb.scaling_equation}")
            if mb.structure_template:
                print(f"           structure: {mb.structure_template}")
    
    if sra_accession:
        print(f"SRA Accession: {sra_accession}")
    
    reverse_map = reverse_base_map(gray_N, gray_A, gray_T, gray_C, gray_G)
    subtract_table = create_base_map()
    inverse_tables = []
    
    inverse_tables = []
    if mode in [1, 2, 3]:
        if metadata_blocks:
            for mb in metadata_blocks:
                formula_func = build_formula_func(mb.scaling_equation)
                inverse_table = build_inverse_quality_table(formula_func, phred_alphabet_max)
                inverse_tables.append(inverse_table)
        
        if not inverse_tables and mode == 1:
            # Mode 1 without metadata - create default inverse table
            inverse_tables.append(np.arange(64, dtype=int))
    
    base_ranges = {
        'A': (1, 63),
        'T': (65, 127),
        'C': (129, 191),
        'G': (193, 255),
        'N': (1, 1)
    }
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    def chunk_generator():
        buffer = data[data_start_byte:]
        chunk_id = 0
        start_seq_idx = 0
        pos = 0
        
        while pos < len(buffer):
            chunk_end = min(pos + chunk_size_bytes, len(buffer))
            
            # Find boundary based on mode
            if mode in [2, 3]:
                # Mode 2: Split on '<' markers
                search_start = max(pos, chunk_end - 1000)
                last_marker = buffer.rfind(b'<', search_start, chunk_end)
                
                if last_marker != -1 and last_marker > pos:
                    seq_end = buffer.find(b'\n', last_marker)
                    if seq_end != -1 and seq_end < len(buffer):
                        chunk_end = seq_end + 1
            else:
                # Mode 0 and 1: Split on '@' headers
                search_start = max(pos, chunk_end - 1000)
                last_header = buffer.rfind(b'\n@', search_start, chunk_end)
                
                if last_header != -1 and last_header > pos:
                    chunk_end = last_header + 1
            
            overlap_start = max(0, pos - 500)
            chunk_binary = buffer[overlap_start:chunk_end]
            
            if chunk_binary:
                actual_start_in_chunk = pos - overlap_start
                if mode == 2:
                    estimated_seqs = chunk_binary[actual_start_in_chunk:].count(b'<')
                else:
                    estimated_seqs = chunk_binary[actual_start_in_chunk:].count(b'\n@')
                
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
    
    with open(output_path, 'wb', buffering=chunk_size_bytes) as outfile:
        with Pool(processes=num_workers) as pool:
            worker_func = partial(
                process_chunk_worker_reconstruction,
                reverse_map=reverse_map,
                subtract_table=subtract_table,
                base_ranges=base_ranges,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                sra_accession=sra_accession,
                mode=mode
            )
            
            for chunk_id, fastq_text, count in pool.imap(worker_func, chunk_generator(), chunksize=1):
                outfile.write(fastq_text)
                
                total_sequences += count
                
                if total_sequences % 100000 == 0:
                    print(f"Reconstructed {total_sequences:,} sequences...")
    
    print(f"\nTotal sequences reconstructed: {total_sequences:,}")
    print(f"Output saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Reconstruct FASTQ files from FASTR format")
    parser.add_argument("input_path", type=str, help="Path to FASTR file")
    parser.add_argument("output_path", type=str, help="Output FASTQ file path")
    
    # Mode selection
    parser.add_argument("--mode", type=int, default=2, choices=[0, 1, 2, 3],
                        help="Conversion mode: 0=headers only, 1=bases only, 2=both, 3=no repeating headers (default: 2)")
    
    # Quality parameters
    parser.add_argument("--phred_offset", type=int, default=33,
                        help="Phred quality offset for output (default: 33)")
    parser.add_argument("--phred_alphabet", type=str, default='phred42',
                        help="Phred alphabet used (phred42/phred63/phred94)")
    
    # Base mapping
    parser.add_argument("--gray_N", type=int, default=1)
    parser.add_argument("--gray_A", type=int, default=63)
    parser.add_argument("--gray_T", type=int, default=127)
    parser.add_argument("--gray_C", type=int, default=191)
    parser.add_argument("--gray_G", type=int, default=255)
    
    # Multiprocessing
    parser.add_argument("--chunk_size_mb", type=int, default=32,
                        help="Chunk size in MB for parallel processing (default: 32)")
    parser.add_argument("--num_workers", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    
    args = parser.parse_args()
    
    if args.phred_alphabet == "phred42":
        phred_alphabet_max = 41
    elif args.phred_alphabet == "phred63":
        phred_alphabet_max = 62
    elif args.phred_alphabet == "phred94":
        phred_alphabet_max = 93
    else:
        phred_alphabet_max = 41
    
    start_time = time.perf_counter()
    
    reconstruct_fastq(
        args.input_path, args.output_path,
        gray_N=args.gray_N, gray_A=args.gray_A,
        gray_T=args.gray_T, gray_C=args.gray_C,
        gray_G=args.gray_G,
        phred_alphabet_max=phred_alphabet_max,
        phred_offset=args.phred_offset,
        chunk_size_mb=args.chunk_size_mb,
        num_workers=args.num_workers,
        mode=args.mode
    )
    
    end_time = time.perf_counter()
    print(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()