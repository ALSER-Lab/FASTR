import os
import tempfile
import numpy as np
import time
import argparse
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from multiprocessing import Pool
from functools import partial
import mmap
import cProfile
import pstats
from io import StringIO
from numba import njit

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
    
    # Find first occurrence of {REPEATING_}
    match = re.search(r'\{REPEATING_\d+\}', structure_template)
    
    if match:
        prefix = structure_template[:match.start()] # Return everything before the first placeholder
        if prefix.endswith(':') or prefix.endswith('/'): # Remove trailing delimiter if present
            prefix = prefix[:-1]
        return prefix
    
    return structure_template


def parse_metadata_header(data: bytes, mode: int) -> Tuple[List[MetadataBlock], int, Optional[str], Optional[int]]:
    """
    Parse metadata headers from the beginning of the file in order to reconstruct it later
    """
    metadata_blocks = []
    sra_accession = None
    phred_alphabet_from_metadata = None
    
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
            if lines and lines[0].startswith('#') and '=' not in lines[0]:
                sra_accession = lines[0][1:]
                print(f"Found SRA accession: {sra_accession}")
                line_idx = 1
            
            # Process metadata lines
            current_mode = None
            current_seq_type = None
            current_structure = None
            current_qual_scale = None
            
            while line_idx < len(lines):
                line = lines[line_idx]
                
                if line.startswith('@'):
                    print(f"Reached sequence headers at line {line_idx}")
                    break
                
                if line.startswith('<'):
                    break
                
                if line.startswith('#MODE='):
                    current_mode = line.split('=', 1)[1].strip()
                    line_idx += 1
                    continue
                
                if line.startswith('#SEQ-TYPE='):
                    current_seq_type = line.split('=', 1)[1].strip()
                    line_idx += 1
                    continue
                
                if line.startswith('#PHRED-ALPHABET='):
                    phred_str = line.split('=', 1)[1].strip()
                    if phred_str.startswith('PHRED_'):
                        try:
                            phred_alphabet_from_metadata = int(phred_str.split('_')[1]) - 1
                            print(f"Found PHRED alphabet: {phred_alphabet_from_metadata}")
                        except:
                            pass
                    line_idx += 1
                    continue
                
                if line.startswith('#STRUCTURE:') or line.startswith('#STRUCTURE='):
                    # Always split on '=' first since that's our actual delimiter
                    if line.startswith('#STRUCTURE='):
                        current_structure = line.split('=', 1)[1].strip()
                    else:
                        current_structure = line.split(':', 1)[1].strip()
                    print(f"Found STRUCTURE metadata: {current_structure}")
                    line_idx += 1
                    continue
                
                if line.startswith('#QUAL_SCALE='):
                    current_qual_scale = line.split('=', 1)[1].strip()
                    print(f"Found equation: {current_qual_scale}")
                    
                    if current_seq_type and current_qual_scale:
                        metadata_blocks.append(MetadataBlock(
                            structure_template=current_structure or "",
                            sequencer_type=current_seq_type,
                            scaling_equation=current_qual_scale,
                            start_index=0,
                            end_index=-1
                        ))
                    
                    line_idx += 1
                    continue
                
                line_idx += 1
        
        # Find actual data start which is the first @ header
        actual_data_start = first_header if first_header > 0 else 0
        return metadata_blocks, actual_data_start, sra_accession, phred_alphabet_from_metadata
    
    # Mode 0, 2, and 3: Header compression enabled
    first_seq_marker = data.find(b'\xff') if mode in [2, 3] else data.find(b'\n@')  # We reserve \xff (255 in hex) for start of sequence indicator
    if first_seq_marker == -1:
        return metadata_blocks, 0, sra_accession, phred_alphabet_from_metadata
    
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
    if lines and lines[0].startswith('#') and '=' not in lines[0]:
        sra_accession = lines[0][1:]
        print(f"Found SRA accession: {sra_accession}")
        line_idx = 1
    
    # Process metadata lines
    current_mode = None
    current_seq_type = None
    current_structure = None
    current_qual_scale = None

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
        
        if line.startswith('#MODE='):
            current_mode = line.split('=', 1)[1].strip()
            line_idx += 1
            continue
        
        if line.startswith('#SEQ-TYPE='):
            current_seq_type = line.split('=', 1)[1].strip()
            line_idx += 1
            continue
        
        if line.startswith('#PHRED-ALPHABET='):
            phred_str = line.split('=', 1)[1].strip()
            if phred_str.startswith('PHRED_'):
                try:
                    phred_alphabet_from_metadata = int(phred_str.split('_')[1]) - 1
                    print(f"Found PHRED alphabet: {phred_alphabet_from_metadata}")
                except:
                    pass
            line_idx += 1
            continue
        
        if line.startswith('#STRUCTURE:') or line.startswith('#STRUCTURE='):
            if line.startswith('#STRUCTURE='):
                current_structure = line.split('=', 1)[1].strip()
            else:
                current_structure = line.split(':', 1)[1].strip()
            print(f"Found STRUCTURE metadata: {current_structure}")
            line_idx += 1
            continue
        
        if line.startswith('#QUAL_SCALE='):
            current_qual_scale = line.split('=', 1)[1].strip()
            print(f"Found equation: {current_qual_scale}")
            line_idx += 1
            continue
        
        line_idx += 1

    if current_seq_type and current_qual_scale:
        metadata_blocks.append(MetadataBlock(
            structure_template=current_structure or "",
            sequencer_type=current_seq_type,
            scaling_equation=current_qual_scale,
            start_index=0,
            end_index=-1
    ))
            
    return metadata_blocks, actual_data_start, sra_accession, phred_alphabet_from_metadata



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

def parse_ont_unique_id(unique_id: str):
    """Parse ONT unique_id which contains key=value pairs"""
    parts = unique_id.strip().split(':')
    
    if parts and '=' not in parts[0]: # First part might not be a kvp but rather a prefix
        prefix = parts[0]
        kv_parts = parts[1:]
    else:
        prefix = ''
        kv_parts = parts
    
    kvs = {}
    for part in kv_parts:
        if '=' in part:
            k, v = part.split('=', 1)
            kvs[k] = v
    
    return prefix, kvs


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

    elif sequencer_type == 'ont':
        result = unique_id
        
        if pair_number > 0:
            result = f"{result}/{pair_number}"
        
        return f"@{result}"
    
    elif sequencer_type == 'ont_sra':
        result = structure.replace('{REPEATING_1}', unique_id)
        
        if pair_number > 0:
            result = f"{result}/{pair_number}"
        
        return f"@{result}"

    else:  # Non-adaptive sequencer types
        delimiter = get_delimiter_for_sequencer(sequencer_type)
        
        # Special handling for PacBio formats that use underscore sub-delimiter
        if sequencer_type in ['pacbio_clr', 'pacbio_subread', 'pacbio_clr_sra']:
            # Split by primary delimiter first
            parts = unique_id.split(delimiter)
            unique_parts = []
            for part in parts:
                # If part contains underscore, split it too
                if '_' in part:
                    unique_parts.extend(part.split('_'))
                else:
                    unique_parts.append(part)
        else:
            unique_parts = unique_id.split(delimiter)
        
        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part)
    
    # Add pair number if present (applies to both adaptive and non-adaptive)
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
    cleaned = re.sub(r'^\s*f\s*\(\s*x\s*\)\s*=\s*', '', formula.strip()).replace('^', '**')

    safe_dict = {
        'ln': np.log, 'log': np.log, 'log10': np.log10, 'exp': np.exp,
        'sqrt': np.sqrt, 'abs': np.abs, 'min': np.minimum, 'max': np.maximum,
        'np': np, '__builtins__': {}
    }

    def formula_func(x):
        local = dict(safe_dict)
        local['x'] = x.astype(np.float32) if isinstance(x, np.ndarray) else np.float32(x)
        with np.errstate(divide='ignore', invalid='ignore'):
            result = eval(cleaned, local)
            if isinstance(result, np.ndarray):
                result = np.nan_to_num(result, nan=1.0, posinf=63.0, neginf=1.0)
            return result

    return formula_func


def create_base_map(gray_N=0, gray_A=3, gray_G=66, gray_C=129, gray_T=192):
    base_table = np.zeros(256, dtype=np.int32)
    base_table[gray_N:gray_A] = gray_N
    base_table[gray_A:gray_G] = gray_A
    base_table[gray_G:gray_C] = gray_G
    base_table[gray_C:gray_T] = gray_C
    base_table[gray_T:255] = gray_T # Never reaches 254 (reserved for indicator of sequence start)
    return base_table


def reverse_base_map(gray_N=0, gray_A=3, gray_G=66, gray_C=129, gray_T=192):
    reverse_map = np.full(256, ord('N'), dtype=np.uint8) # Default to 'N'
    reverse_map[gray_N:gray_A] = ord('N')
    reverse_map[gray_A:gray_G] = ord('A')
    reverse_map[gray_G:gray_C] = ord('G')
    reverse_map[gray_C:gray_T] = ord('C')
    reverse_map[gray_T:255] = ord('T') # Never reaches 254 (reserved for indicator of sequence start)
    return reverse_map

@njit
def reverse_scaling_to_quality(binary_values, subtract_table, inverse_table, max_phred):
    y = binary_values - subtract_table[binary_values]
    
    # Clip between 0 and 63
    y_clipped = np.empty_like(y)
    for i in range(y.shape[0]):
        val = y[i]
        if val < 0:
            val = 0
        elif val > 63:
            val = 63
        y_clipped[i] = val
    
    # Map through inverse table
    result = np.empty_like(y_clipped)
    for i in range(y_clipped.shape[0]):
        val = inverse_table[y_clipped[i]]
        if val < 0:
            val = 0
        elif val > max_phred:
            val = max_phred
        result[i] = val
    
    return result

@njit
def build_inverse_quality_table(scaled_int_for_q, q_possible, max_range):
    inverse_table = np.full(max_range + 1, -1, dtype=np.int32)
    
    for i in range(q_possible.shape[0]):
        s = scaled_int_for_q[i]
        if 0 <= s <= max_range:
            inverse_table[s] = int(q_possible[i])
    
    # Fill gaps
    valid_indices = []
    for i in range(max_range + 1):
        if inverse_table[i] != -1:
            valid_indices.append(i)
    
    if len(valid_indices) > 0:
        # Fill beginning
        for i in range(valid_indices[0]):
            inverse_table[i] = inverse_table[valid_indices[0]]
        
        # Fill gaps
        for i in range(len(valid_indices) - 1):
            start_idx = valid_indices[i]
            end_idx = valid_indices[i + 1]
            start_val = inverse_table[start_idx]
            end_val = inverse_table[end_idx]
            gap = end_idx - start_idx
            for j in range(1, gap):
                inverse_table[start_idx + j] = int(round(start_val + (end_val - start_val) * j / gap))
        
        # Fill end
        for i in range(valid_indices[-1] + 1, max_range + 1):
            inverse_table[i] = inverse_table[valid_indices[-1]]
    else:
        for i in range(max_range + 1):
            inverse_table[i] = 0
    
    return inverse_table

def quality_to_ascii(quality_scores: np.ndarray, phred_offset: int = 33) -> bytes:
    """Convert numeric quality scores to ASCII string"""
    return bytes((quality_scores + phred_offset).astype(np.uint8))


def process_chunk_worker_reconstruction(chunk_data, mmap_path, reverse_map, subtract_table, 
                                       base_ranges, metadata_blocks, inverse_tables, 
                                       phred_alphabet_max, phred_offset, sra_accession,
                                       mode, headers_file_path=None):
    """
    Worker function that processes a single chunk of binary sequence data in parallel.
    For mode 3, loads headers from file instead of receiving them via pickle.
    """
    try:
        chunk_id, abs_start, abs_end, start_seq_idx = chunk_data
        
        # For mode 3, memory-map headers file in worker (avoiding pickle overhead)
        headers_mmap = None
        headers_offsets = None
        if mode == 3 and headers_file_path:
            with open(headers_file_path, 'rb') as hf:
                headers_mmap = mmap.mmap(hf.fileno(), 0, access=mmap.ACCESS_READ)
                # Build line offset index for this chunk's range only
                # This is MUCH faster than splitting the entire file!!
                headers_offsets = []
                pos = 0
                line_idx = 0
                # We only need offsets for sequences in this chunk
                target_start = start_seq_idx
                target_end = start_seq_idx + 200000  # Estimate, will stop when done
                
                while pos < len(headers_mmap) and line_idx < target_end:
                    if line_idx >= target_start:
                        headers_offsets.append(pos)
                    next_newline = headers_mmap.find(b'\n', pos)
                    if next_newline == -1:
                        if line_idx >= target_start:
                            headers_offsets.append(pos)
                        break
                    pos = next_newline + 1
                    line_idx += 1
        
        # Open shared mmap in worker
        with open(mmap_path, 'rb') as f:
            data = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            chunk_binary = data[abs_start:abs_end]
        
        temp = tempfile.NamedTemporaryFile(delete=False, mode='wb', buffering=8*1024*1024)
        sequence_count = start_seq_idx
        
        # Determine which metadata block to use for this chunk
        current_metadata_idx = 0
        if metadata_blocks:
            current_metadata = metadata_blocks[0]
            if inverse_tables:
                current_inverse_table = inverse_tables[0]
            else:
                current_inverse_table = None
        else:
            current_metadata = None
            current_inverse_table = None
        
        # Mode 1: Only bases compressed (original headers intact)
        if mode == 1:
            cursor = 0
            while cursor < len(chunk_binary):
                # Find the start of a header
                header_start = chunk_binary.find(b'@', cursor)
                if header_start == -1: break
                
                # Original header
                header_end = chunk_binary.find(b'\n', header_start)
                if header_end == -1: break
                header = chunk_binary[header_start:header_end].decode('utf-8', errors='ignore')
                
                # Look for binary marker \xff after the header
                marker_pos = chunk_binary.find(b'\xff', header_end)
                if marker_pos != -1:
                    # Binary sequence
                    seq_start = marker_pos + 1
                    
                    # Find where this record ends (the next @ or end of chunk)
                    next_record = chunk_binary.find(b'\n@', seq_start)
                    seq_end = next_record if next_record != -1 else len(chunk_binary)
                    
                    seq_data = chunk_binary[seq_start:seq_end].rstrip(b'\r\n') # Skip 255 marker
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
                        quality_string = 'I' * len(bases) # If you see all I's (even when there shouldn't be) it means we had to fallback
                    
                    temp.write(f"{header}\n{bases}\n+\n{quality_string}\n".encode('ascii'))
                    sequence_count += 1
                    cursor = seq_end
                else:
                    cursor = header_end + 1

        elif mode == 0:
            cursor = 0
            while cursor < len(chunk_binary):
                # Find start of record
                start_idx = chunk_binary.find(b'@', cursor)
                if start_idx == -1: break
                
                header_end = chunk_binary.find(b'\n', start_idx)
                if header_end == -1: break
                
                header_content = chunk_binary[start_idx+1:header_end].decode('utf-8', errors='ignore').strip()
                
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
                
                plus_idx = chunk_binary.find(b'\n+\n', header_end)
                if plus_idx == -1: break
                
                bases_data = chunk_binary[header_end:plus_idx].strip()
                bases = bases_data.decode('ascii', errors='ignore')
                
                qual_start = plus_idx + 3
                qual_data = chunk_binary[qual_start : qual_start + len(bases_data)]
                quality_string = qual_data.decode('ascii', errors='ignore')
                
                temp.write(f"{header}\n{bases}\n+\n{quality_string}\n".encode('ascii'))
                sequence_count += 1
                cursor = qual_start + len(bases_data)

        if mode == 3:
            # Find all markers
            chunk_array = np.frombuffer(chunk_binary, dtype=np.uint8)
            marker_positions = np.where(chunk_array == 255)[0].tolist()
            
            for idx, marker_pos in enumerate(marker_positions):
                # Binary sequence
                seq_start = marker_pos + 1 # Skip 255 marker
                
                # Sequence ends at the next marker
                if idx < len(marker_positions) - 1:
                    seq_end = marker_positions[idx+1]
                else:
                    seq_end = len(chunk_binary)
                
                seq_data = chunk_binary[seq_start:seq_end].rstrip(b'\r\n')
                
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
                
                if headers_mmap and headers_offsets: # Get using mmap and offset index
                    offset_idx = sequence_count - start_seq_idx
                    if offset_idx < len(headers_offsets):
                        start_offset = headers_offsets[offset_idx]
                        if offset_idx + 1 < len(headers_offsets):
                            end_offset = headers_offsets[offset_idx + 1] - 1  # -1 for newline
                        else:
                            # Last header in our range
                            end_offset = headers_mmap.find(b'\n', start_offset)
                            if end_offset == -1:
                                end_offset = len(headers_mmap)
                        
                        header = headers_mmap[start_offset:end_offset].decode('utf-8', errors='ignore')
                        if not header.startswith('@'):
                            header = '@' + header
                    else:
                        # Fallback
                        header = f"@seq{sequence_count}"
                else:
                    # Fallback
                    if current_metadata and current_metadata.structure_template: 
                        header_prefix = find_structure_prefix(current_metadata.structure_template)
                        header = f"@{header_prefix}.{sequence_count}" if header_prefix else f"@seq{sequence_count}"
                    elif sra_accession:
                        header = f"@{sra_accession}.{sequence_count}"
                    else:
                        header = f"@seq{sequence_count}"
                
                # Reconstruct quality
                if current_inverse_table is not None:
                    quality_scores = reverse_scaling_to_quality(
                        seq_array, subtract_table,
                        current_inverse_table, max_phred=phred_alphabet_max
                    )
                    quality_string = bytes((quality_scores + phred_offset).astype(np.uint8)).decode('ascii')
                else:
                    quality_string = 'I' * len(bases)
                
                temp.write(f"{header}\n".encode('ascii'))
                temp.write(bases.encode('ascii'))
                temp.write(b"\n+\n")
                temp.write(quality_string.encode('ascii'))
                temp.write(b"\n")
                sequence_count += 1

        # Mode 2: Headers compressed, bases as binary
        else:
            sequences_in_chunk = 0
            
            # Find all 255 markers in this chunk first
            chunk_array = np.frombuffer(chunk_binary, dtype=np.uint8)
            marker_positions = np.where(chunk_array == 255)[0].tolist()
            
            for marker_idx, marker_pos in enumerate(marker_positions):
                # Find the header before this marker
                line_start = chunk_binary.rfind(b'\n@', max(0, marker_pos - 500), marker_pos)
                
                if line_start == -1:
                    # Check if file/chunk starts with @
                    if chunk_binary[:1] == b'@':
                        line_start = -1  # Will become 0 after +1
                    else: # Skip if no header found
                        continue
                
                header_start = line_start + 2  # Skip '\n@'
                if line_start == -1:
                    header_start = 1  # Skip just '@' at start of file
                
                # Find the end of the header line (the \n before 255 marker)
                header_end = chunk_binary.rfind(b'\n', header_start, marker_pos)
                if header_end == -1 or header_end < header_start:
                    header_end = marker_pos
                
                header_content = chunk_binary[header_start:header_end].decode('utf-8', errors='ignore').strip()
                
                # Parse pair number if present
                unique_id = None
                pair_number = 0
                
                if '/' in header_content:
                    parts = header_content.rsplit('/', 1)
                    # Only treat as pair number if the last part is a single digit (1 or 2)
                    if parts[1].isdigit() and len(parts[1]) == 1 and int(parts[1]) in [1, 2]:
                        unique_id = parts[0]
                        pair_number = int(parts[1])
                    else:
                        # Not a pair number, use full header_content as unique_id
                        unique_id = header_content
                else:
                    unique_id = header_content
                
                # Reconstruct the full header
                if unique_id and current_metadata and current_metadata.structure_template:
                    header = reconstruct_header_from_structure(
                        current_metadata.structure_template,
                        unique_id,
                        current_metadata.sequencer_type,
                        pair_number
                    )
                else:
                    if sra_accession:
                        header = f"@{sra_accession}.{sequence_count}"
                    else:
                        header = f"@seq{sequence_count}"
                
                # New debug: cannot search for \n as an indicator because binary data represents '\n' as 10!!
                # This was causing premature stopping
                seq_start = marker_pos + 1  # Start after 255 marker
                
                if marker_idx < len(marker_positions) - 1:
                    # Not the last marker - sequence ends before next marker
                    next_marker = marker_positions[marker_idx + 1]
                    # Find the newline before the next marker (which precedes the next header)
                    seq_end_search = chunk_binary.rfind(b'\n@', seq_start, next_marker)
                    if seq_end_search != -1:
                        seq_end = seq_end_search
                    else:
                        # No header found, use position right before next marker
                        seq_end = next_marker
                else:
                    # Last marker - sequence goes to end of chunk or next header
                    next_header = chunk_binary.find(b'\n@', seq_start)
                    if next_header != -1:
                        seq_end = next_header
                    else:
                        # No next header, use end of chunk
                        seq_end = len(chunk_binary)
                
                seq_data = chunk_binary[seq_start:seq_end]
                
                # Remove trailing newline if present
                if seq_data and seq_data[-1:] == b'\n':
                    seq_data = seq_data[:-1]
                
                if len(seq_data) == 0: # Skip if sequence data is empty
                    continue
                
                # Check for metadata block changes
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_inverse_table = inverse_tables[current_metadata_idx + 1]
                        current_metadata_idx += 1
                
                seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                bases_array = reverse_map[seq_array]
                bases = bases_array.tobytes().decode('ascii')
                
                # Reconstruct quality
                if current_inverse_table is not None:
                    quality_scores = reverse_scaling_to_quality(
                        seq_array, subtract_table,
                        current_inverse_table, max_phred=phred_alphabet_max
                    )
                    quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
                else:
                    quality_string = 'I' * len(bases)

                temp.write(header.encode('ascii'))
                temp.write(b'\n')
                temp.write(bases.encode('ascii'))
                temp.write(b'\n+\n')
                temp.write(quality_string.encode('ascii'))
                temp.write(b'\n')
                sequence_count += 1
                sequences_in_chunk += 1
            
        
        # print(f"Worker completed chunk {chunk_id}, processed {sequence_count - start_seq_idx} sequences")
        temp.close()
        return temp.name, sequence_count - start_seq_idx
    
    except Exception as e:
        print(f"ERROR in worker processing chunk {chunk_id}: {e}")
        import traceback
        traceback.print_exc()
        raise


def reconstruct_fastq(input_path: str, output_path: str, 
                     gray_N: int = 0, gray_A: int = 3, gray_G: int = 66,
                     gray_C: int = 129, gray_T: int = 192, phred_alphabet_max: int = None,
                     phred_offset: int = 33, chunk_size_mb: int = 8, num_workers: int = 4,
                     mode: int = 2, mode3_headers_file: str = None):
    """
    Reconstruct FASTQ file from FASTR using parallel processing.
    """
    print(f"Reading FASTR: {input_path}")
    print(f"Reconstruction mode: {mode}")
    
    # For mode 3, validate headers file exists but don't load it..
    # Workers will load it themselves to avoid pickle overhead
    if mode == 3 and mode3_headers_file:
        if not os.path.exists(mode3_headers_file):
            raise FileNotFoundError(f"Headers file not found: {mode3_headers_file}")
        print(f"Headers file: {mode3_headers_file} (will be loaded by workers)")
        # Quick check of header count
        with open(mode3_headers_file, 'rb') as hf:
            header_count = sum(1 for _ in hf)
        print(f"Headers file contains {header_count:,} headers")
    elif mode == 3 and not mode3_headers_file:
        print("WARNING: Mode 3 requires --headers_file argument")
        print("Proceeding without headers - will use fallback header generation (currently: '@seq')")
    
    with open(input_path, 'rb') as f:
        data = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    file_size = len(data)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    print(f"Using {num_workers} worker processes for parallel reconstruction")
    
    metadata_blocks, data_start_byte, sra_accession, phred_from_metadata = parse_metadata_header(data, mode)
    
    if phred_alphabet_max is None:
        if phred_from_metadata is not None:
            phred_alphabet_max = phred_from_metadata
            print(f"Using PHRED alphabet from metadata: {phred_alphabet_max}")
        else:
            phred_alphabet_max = 41
            print(f"No PHRED alphabet found in metadata, using default: {phred_alphabet_max}")
    else:
        print(f"Using user-specified PHRED alphabet (overriding metadata): {phred_alphabet_max}")
    
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
    
    reverse_map = reverse_base_map(gray_N, gray_A, gray_G, gray_C, gray_T)
    subtract_table = create_base_map(gray_N, gray_A, gray_G, gray_C, gray_T)
    
    inverse_tables = []
    if mode in [1, 2, 3]:
        if metadata_blocks:
            for mb in metadata_blocks:
                formula_func = build_formula_func(mb.scaling_equation)
                q_possible = np.arange(phred_alphabet_max + 1, dtype=np.float32)
                scaled_int_for_q = formula_func(q_possible).astype(np.int32)
                inverse_table = build_inverse_quality_table(
                    scaled_int_for_q, 
                    q_possible.astype(np.int32), 
                    63  # max_range should be 63 as used elsewhere
                )
                inverse_tables.append(inverse_table)
        
        if not inverse_tables and mode == 1:
            # Mode 1 without metadata, so we'll create default inverse table
            inverse_tables.append(np.arange(64, dtype=int))
    
    base_ranges = {
        'A': (3, 65),
        'C': (66, 128),
        'G': (129, 191),
        'T': (192, 254),
        'N': (0, 2)
    }
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    def chunk_generator():
        buffer_start = data_start_byte
        buffer_len = len(data) - data_start_byte
        chunk_id = 0
        start_seq_idx = 0
        pos = 0
        
        while pos < buffer_len:
            chunk_end = min(pos + chunk_size_bytes, buffer_len)
            
            # Find boundary based on mode
            if mode in [2, 3]:
                # Mode 2/3: Split on 255 markers
                search_start = max(pos, chunk_end - 1000)
                last_marker = data[buffer_start + search_start:buffer_start + chunk_end].rfind(b'\xff')
                
                if last_marker != -1 and (search_start + last_marker) > pos:
                    seq_end_offset = data[buffer_start + search_start + last_marker:buffer_start + chunk_end].find(b'\n')
                    if seq_end_offset != -1:
                        chunk_end = search_start + last_marker + seq_end_offset + 1
            else:
                # Mode 0 and 1: Split on '@' headers
                search_start = max(pos, chunk_end - 1000)
                last_header = data[buffer_start + search_start:buffer_start + chunk_end].rfind(b'\n@')
                
                if last_header != -1 and (search_start + last_header) > pos:
                    chunk_end = search_start + last_header + 1
            
            overlap_start = max(0, pos - 500)
            abs_start = buffer_start + overlap_start
            abs_end = buffer_start + chunk_end
            
            if abs_end > abs_start:
                actual_start_in_chunk = pos - overlap_start
                chunk_view = data[abs_start:abs_end]
                if mode in [2, 3]:
                    # Count 255 markers for mode 2/3
                    estimated_seqs = chunk_view[actual_start_in_chunk:].count(b'\xff')
                else:
                    estimated_seqs = chunk_view[actual_start_in_chunk:].count(b'\n@')
                
                yield (chunk_id, abs_start, abs_end, start_seq_idx)
                
                start_seq_idx += estimated_seqs
                chunk_id += 1
                if chunk_id % 10 == 0:
                    print(f"Read {chunk_id} chunks...")
            
            pos = chunk_end
            if pos >= buffer_len:
                break
    
    print("Reconstructing FASTQ with parallel processing...")
    total_sequences = 0
    
    with open(output_path, 'wb', buffering=chunk_size_bytes) as outfile:
        with Pool(processes=num_workers) as pool:
            worker_func = partial( # Redid argument order as it was causing some np arrays to be pickled instead of metadata objects
                process_chunk_worker_reconstruction,
                mmap_path=input_path,        
                reverse_map=reverse_map,       
                subtract_table=subtract_table,  
                base_ranges=base_ranges,     
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                sra_accession=sra_accession,
                mode=mode,
                headers_file_path=mode3_headers_file  # Pass file path instead of list
            )
            
            temp_files = []
            for temp_path, count in pool.imap(worker_func, chunk_generator(), chunksize=1):
                temp_files.append(temp_path)
                total_sequences += count
                if total_sequences % 100000 == 0:
                    print(f"Reconstructed {total_sequences:,} sequences...")
        
        # Merge temp files into the final output and clean up
        for temp_path in temp_files:
            with open(temp_path, 'rb') as tf:
                outfile.write(tf.read())
            os.remove(temp_path)

    
    print(f"\nTotal sequences reconstructed: {total_sequences:,}")
    print(f"Output saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Reconstruct FASTQ files from FASTR format")
    parser.add_argument("input_path", type=str, help="Path to FASTR file")
    parser.add_argument("output_path", type=str, help="Output FASTQ file path")
    
    # Mode selection
    parser.add_argument("--mode", type=int, default=2, choices=[0, 1, 2, 3],
                        help="Conversion mode: 0=headers only, 1=bases only, 2=both, 3=no repeating headers (default: 2)")
    parser.add_argument("--headers_file", type=str, default=None,
                        help="Path to headers file for mode 3 reconstruction")
    
    
    # Quality parameters
    parser.add_argument("--phred_offset", type=int, default=33,
                        help="Phred quality offset for output (default: 33)")
    parser.add_argument("--phred_alphabet", type=str, default=None,
                        help="Phred alphabet override (phred42/phred63/phred94), defaults to metadata value")
    parser.add_argument("--profile", action="store_true",
                        help="Enable cProfile profiling")
    
    # Base mapping
    parser.add_argument("--gray_N", type=int, default=0)
    parser.add_argument("--gray_A", type=int, default=3)
    parser.add_argument("--gray_C", type=int, default=129)
    parser.add_argument("--gray_G", type=int, default=66)
    parser.add_argument("--gray_T", type=int, default=192)
    
    # Multiprocessing
    parser.add_argument("--chunk_size_mb", type=int, default=8,
                        help="Chunk size in MB for parallel processing (default: 8)")
    parser.add_argument("--num_workers", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    
    args = parser.parse_args()
    
    phred_alphabet_max = None
    if args.phred_alphabet:
        if args.phred_alphabet == "phred42":
            phred_alphabet_max = 41
        elif args.phred_alphabet == "phred63":
            phred_alphabet_max = 62
        elif args.phred_alphabet == "phred94":
            phred_alphabet_max = 93
    
    start_time = time.perf_counter()

    if args.profile:
        profiler = cProfile.Profile()
        profiler.enable()
    
    reconstruct_fastq(
        args.input_path, args.output_path,
        gray_N=args.gray_N, gray_A=args.gray_A,
        gray_T=args.gray_T, gray_C=args.gray_C,
        gray_G=args.gray_G,
        phred_alphabet_max=phred_alphabet_max,
        phred_offset=args.phred_offset,
        chunk_size_mb=args.chunk_size_mb,
        num_workers=args.num_workers,
        mode=args.mode,
        mode3_headers_file=args.headers_file
    )
    if args.profile:
        profiler.disable()
        s = StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
        ps.print_stats(20) # Limit output so it doesnt take up all space in terminal
        print(s.getvalue())
    
    end_time = time.perf_counter()
    print(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()