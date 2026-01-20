import io
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
import cProfile
import pstats
from io import StringIO
from numba import njit
import traceback
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

__version__ = "1.0.0"

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
    length_flag = False
    second_head_flag = False
    phred_alphabet_from_metadata = None
    detected_mode = mode
    safe_mode_flag = False
    
    if mode == 1:
        # Mode 1: Only bases compressed (no header compression)
        # Parse metadata for quality reconstruction, then find where actual sequences start

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
            
            logger.info(f"Parsed {len(lines)} header lines")
            for i, line in enumerate(lines[:5]): 
                logger.info(f"  Line {i}: {line[:80]}")
            
            line_idx = 0
            
            # Check for SRA accession
            if lines and lines[0].startswith('#') and '=' not in lines[0]:
                sra_accession = lines[0][1:]
                logger.info(f"Found SRA accession: {sra_accession}")
                line_idx = 1
            
            # Process metadata lines
            current_seq_type = None
            current_structure = None
            current_qual_scale = None
            
            while line_idx < len(lines):
                line = lines[line_idx]
                
                if line.startswith('@'):
                    logger.info(f"Reached sequence headers at line {line_idx}")
                    break
                
                if line.startswith('<'):
                    break
                
                if line.startswith('#MODE='):
                    detected_mode = int(line.split('=', 1)[1].strip())
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
                            logger.info(f"Found PHRED alphabet: {phred_alphabet_from_metadata}")
                        except:
                            pass
                    line_idx += 1
                    continue

                if line.startswith('#LENGTH='):
                    length_str = line.split('=', 1)[1].strip()
                    if length_str.lower() == 'y':
                        length_flag = True
                        logger.info(f"Found LENGTH flag: {length_flag}")
                    line_idx += 1
                    continue

                if line.startswith('#SECOND_HEAD='):
                    second_head_str = line.split('=', 1)[1].strip()
                    if second_head_str.lower() == 'y':
                        second_head_flag = True
                        logger.info(f"Found SECOND_HEAD flag: {second_head_flag}")
                    line_idx += 1
                    continue
                
                if line.startswith('#SAFE_MODE='):
                    safe_mode_str = line.split('=', 1)[1].strip()
                    if safe_mode_str.lower() == 'y':
                        safe_mode_flag = True
                        logger.info(f"Found SAFE_MODE flag: {second_head_flag}")
                    line_idx += 1
                    continue

                if line.startswith('#STRUCTURE:') or line.startswith('#STRUCTURE='):
                    # Always split on '=' first since that's our actual delimiter
                    if line.startswith('#STRUCTURE='):
                        current_structure = line.split('=', 1)[1].strip()
                    else:
                        current_structure = line.split(':', 1)[1].strip()
                    logger.info(f"Found STRUCTURE metadata: {current_structure}")
                    line_idx += 1
                    continue
                
                if line.startswith('#QUAL_SCALE='):
                    current_qual_scale = line.split('=', 1)[1].strip()
                    logger.info(f"Found equation: {current_qual_scale}")
                    
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
        return metadata_blocks, actual_data_start, sra_accession, phred_alphabet_from_metadata, detected_mode, length_flag, second_head_flag, safe_mode_flag
    
    # Mode 0, 2, and 3: Header compression enabled
    # We reserve \xff (255 in hex) for start of sequence indicator, only for mode 3 (given it doesn't have an '@' indicator)
    first_seq_marker = data.find(b'\xff') if mode == 3 else data.find(b'\n@') 
    if first_seq_marker == -1:
        return metadata_blocks, 0, sra_accession, phred_alphabet_from_metadata, detected_mode, length_flag, second_head_flag, safe_mode_flag
    
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
    
    logger.info(f"Parsed {len(lines)} header lines")
    for i, line in enumerate(lines[:5]):  # Show first 5 lines
        logger.info(f"  Line {i}: {line[:80]}")
    
    line_idx = 0
    
    # Check for SRA accession
    if lines and lines[0].startswith('#') and '=' not in lines[0]:
        sra_accession = lines[0][1:]
        logger.info(f"Found SRA accession: {sra_accession}")
        line_idx = 1
    
    # Process metadata lines
    current_seq_type = None
    current_structure = None
    current_qual_scale = None

    while line_idx < len(lines):
        line = lines[line_idx]
        
        if line.startswith('<'):
            break
        
        if line.startswith('#MODE='):
            detected_mode = int(line.split('=', 1)[1].strip())
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
                    logger.info(f"Found PHRED alphabet: {phred_alphabet_from_metadata}")
                except:
                    pass
            line_idx += 1
            continue
        
        if line.startswith('#LENGTH='):
            length_str = line.split('=', 1)[1].strip()
            if length_str.lower() == 'y':
                length_flag = True
                logger.info(f"Found LENGTH flag: {length_flag}")
            line_idx += 1
            continue

        if line.startswith('#SECOND_HEAD='):
            second_head_str = line.split('=', 1)[1].strip()
            if second_head_str.lower() == 'y':
                second_head_flag = True
                logger.info(f"Found SECOND_HEAD flag: {second_head_flag}")
            line_idx += 1
            continue

        if line.startswith('#SAFE_MODE='):
            safe_mode_str = line.split('=', 1)[1].strip()
            if safe_mode_str.lower() == 'y':
                safe_mode_flag = True
                logger.info(f"Found SAFE_MODE flag: {second_head_flag}")
            line_idx += 1
            continue
        
        if line.startswith('#STRUCTURE:') or line.startswith('#STRUCTURE='):
            if line.startswith('#STRUCTURE='):
                current_structure = line.split('=', 1)[1].strip()
            else:
                current_structure = line.split(':', 1)[1].strip()
            logger.info(f"Found STRUCTURE metadata: {current_structure}")
            line_idx += 1
            continue
        
        if line.startswith('#QUAL_SCALE='):
            current_qual_scale = line.split('=', 1)[1].strip()
            logger.info(f"Found equation: {current_qual_scale}")
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
            
    return metadata_blocks, actual_data_start, sra_accession, phred_alphabet_from_metadata, detected_mode, length_flag, second_head_flag, safe_mode_flag

def get_delimiter_for_sequencer(sequencer_type: str) -> str:
    """Get the delimiter character used by each sequencer type"""
    if sequencer_type == 'illumina' or sequencer_type == 'old_illumina':
        return ':'
    elif sequencer_type.startswith('pacbio') and not sequencer_type.endswith('_sra'):
        return '/'
    elif sequencer_type == 'ont':
        return ':'
    elif sequencer_type == 'sra' or sequencer_type.endswith('_sra'):
        return ' ' 
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
    
    elif sequencer_type == 'illumina_sra':
        unique_parts_by_space = unique_id.split(' ')
        
        if len(unique_parts_by_space) < 2:
            result = structure.replace('{REPEATING_1}', unique_id)
        else:
            spot = unique_parts_by_space[0]
            illumina_coords = unique_parts_by_space[1]  
            coord_parts = illumina_coords.split(':')
            
            result = structure.replace('{REPEATING_1}', spot)
            
            for i, coord_part in enumerate(coord_parts):
                placeholder = f"{{REPEATING_{i+2}}}" 
                result = result.replace(placeholder, coord_part, 1)
            
            if len(unique_parts_by_space) > 2:
                extra_fields = ' '.join(unique_parts_by_space[2:])
                result = result.replace('{REPEATING_5}', extra_fields, 1)
    
    elif sequencer_type in ['pacbio_hifi_sra', 'pacbio_clr_sra']:
        parts_by_space = unique_id.split(' ')
        
        if len(parts_by_space) < 2:
            # Fallback
            result = structure
            for i, part in enumerate(parts_by_space, 1):
                result = result.replace(f"{{REPEATING_{i}}}", part, 1)
        else:
            spot = parts_by_space[0]
            pacbio_path = parts_by_space[1]  
            result = structure.replace('{REPEATING_1}', spot, 1)
            path_parts = pacbio_path.split('/')
            
            if sequencer_type == 'pacbio_clr_sra' and len(path_parts) == 2:
                hole = path_parts[0]
                coords = path_parts[1]
                
                result = result.replace('{REPEATING_2}', hole, 1)
                
                if '_' in coords:
                    coord_parts = coords.split('_')
                    result = result.replace('{REPEATING_3}', coord_parts[0], 1)
                    if len(coord_parts) > 1:
                        result = result.replace('{REPEATING_4}', coord_parts[1], 1)
                else:
                    result = result.replace('{REPEATING_3}', coords, 1)
            
            elif sequencer_type == 'pacbio_hifi_sra':
                if len(path_parts) >= 1:
                    result = result.replace('{REPEATING_2}', path_parts[0], 1)
            
            if len(parts_by_space) > 2:
                extra = ' '.join(parts_by_space[2:])
                result = result.replace('{REPEATING_5}', extra, 1)  
                result = result.replace('{REPEATING_3}', extra, 1)  
    
    elif sequencer_type == 'sra':
        unique_parts = unique_id.split(' ')
        result = structure
        for i, part in enumerate(unique_parts, 1):
            placeholder = f"{{REPEATING_{i}}}"
            result = result.replace(placeholder, part, 1)

    else:  # Non-adaptive, non-SRA sequencer types
        delimiter = get_delimiter_for_sequencer(sequencer_type)
        
        # Special handling for PacBio formats that use underscore sub-delimiter
        if sequencer_type in ['pacbio_clr', 'pacbio_subread']:
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
    
    # Add pair number if present (applies to all formats)
    if pair_number > 0:
        result = f"{result}/{pair_number}"
    
    return f"@{result}"


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

def build_header_index(headers_file_path: str) -> str:
    logger.info(f"Building header index for {headers_file_path}...")
    
    line_count = 0
    with open(headers_file_path, 'rb') as hf:
        for _ in hf:
            line_count += 1
    
    logger.info(f"Found {line_count:,} headers, creating index...")
    
    temp_dir = tempfile.gettempdir()
    mmap_path = os.path.join(temp_dir, f'header_index_{os.getpid()}.mmap')
    
    offsets_mmap = np.memmap(mmap_path, dtype=np.int64, mode='w+', shape=(line_count,))
    
    with open(headers_file_path, 'rb') as hf:
        offset = 0
        idx = 0
        while True:
            line = hf.readline()
            if not line:
                break
            offsets_mmap[idx] = offset
            offset = hf.tell()
            idx += 1
    
    offsets_mmap.flush()
    logger.info(f"Header index built: {line_count:,} headers")
    
    return offsets_mmap, mmap_path


def get_header_by_index(headers_file_path: str, index_array: np.ndarray, seq_idx: int) -> str:
    if seq_idx >= len(index_array):
        return None
    
    offset = int(index_array[seq_idx])  
    with open(headers_file_path, 'rb') as hf:
        hf.seek(offset)
        header_line = hf.readline()
        return header_line.decode('utf-8', errors='ignore').strip()



def process_chunk_worker_reconstruction(chunk_data, mmap_path, reverse_map, subtract_table, 
                                       metadata_blocks, inverse_tables, 
                                       phred_alphabet_max, phred_offset, sra_accession,
                                       mode, length_flag=False, headers_file_path=None, headers_mmap_info=None,
                                       second_head_flag=False, safe_mode_flag=False, data_start_byte=None):
    """
    Worker function that processes a single chunk of binary sequence data.
    """
    try:
        chunk_id, abs_start, abs_end, start_seq_idx = chunk_data
        
        if mode == 3 and headers_mmap_info:
            if not hasattr(process_chunk_worker_reconstruction, '_header_index'):
                mmap_path_hdr, shape, dtype = headers_mmap_info
                process_chunk_worker_reconstruction._header_index = np.memmap(
                    mmap_path_hdr, dtype=dtype, mode='r', shape=shape
                )
                logger.debug(f"Worker loaded header index: {shape[0]:,} offsets")
            header_index = process_chunk_worker_reconstruction._header_index
        else:
            header_index = None
        
        MAX_CHUNK_EXTENSION = 10 * 1024 * 1024
        
        with open(mmap_path, 'rb') as f:
            f.seek(abs_start)
            chunk_read_size = (abs_end - abs_start) + MAX_CHUNK_EXTENSION
            data = f.read(chunk_read_size)
        
        output_buffer = io.BytesIO()
        sequence_count = start_seq_idx
        local_count = 0
        
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
            sequences_in_chunk = 0
            cursor = 0
            chunk_size = abs_end - abs_start
            
            while cursor < chunk_size:
                if safe_mode_flag:
                    # Safe mode, format is: \n@header(255)\nbases\n
                    # Find next @ with validation
                    if cursor == 0 and data[0:1] == b'@':
                        # First sequence in chunk starts at position 0
                        header_start_rel = 0
                    else:
                        # Always look for \n@ for subsequent sequences
                        header_start_rel = data.find(b'\n@', cursor)
                        if header_start_rel != -1:
                            header_start_rel += 1
                    
                    if header_start_rel < 0 or header_start_rel >= chunk_size:
                        break
                    
                    candidate = header_start_rel # Then validate by finding \xff after this
                    
                    while candidate >= 0:
                        next_xff_rel = data.find(b'\xff', candidate, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                        next_at_rel = data.find(b'\n@', candidate + 1, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                        # If another \n@ comes before \xff, current candidate is invalid
                        if next_at_rel != -1 and (next_xff_rel == -1 or next_at_rel < next_xff_rel):
                            candidate = next_at_rel + 1
                            if candidate >= chunk_size + MAX_CHUNK_EXTENSION:
                                candidate = -1
                                break
                            continue
                        
                        if next_xff_rel != -1:
                            xff_pos = next_xff_rel
                            break
                        
                        candidate = -1
                        break
                    
                    if candidate == -1:
                        break
                    
                    header_start = candidate
                    
                    header = data[header_start:xff_pos].decode('utf-8', errors='replace').strip()
                    
                    if xff_pos + 1 < len(data) and data[xff_pos+1:xff_pos+2] == b'\n':
                        seq_start_rel = xff_pos + 2
                    else:
                        seq_start_rel = xff_pos + 1
                    
                    candidate_end_rel = seq_start_rel
                    seq_end_rel = -1
                    
                    while True:
                        next_at_rel = data.find(b'\n@', candidate_end_rel, min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION))
                        
                        if next_at_rel != -1:
                            next_xff_rel = data.find(b'\xff', next_at_rel + 1, min(len(data), next_at_rel + 1000))
                            following_at_rel = data.find(b'\n@', next_at_rel + 2, min(len(data), next_at_rel + 1000))
                            
                            if following_at_rel != -1 and (next_xff_rel == -1 or following_at_rel < next_xff_rel):
                                candidate_end_rel = following_at_rel + 1
                                if candidate_end_rel > seq_start_rel + MAX_CHUNK_EXTENSION:
                                    break
                                continue
                            
                            seq_end_rel = next_at_rel
                            break
                        else:
                            seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                            while seq_end_rel > seq_start_rel and data[seq_end_rel-1:seq_end_rel] in (b'\n', b'\r', b' '):
                                seq_end_rel -= 1
                            break
                    
                    if seq_end_rel == -1 or seq_end_rel <= seq_start_rel:
                        cursor = len(data)
                        continue
                    
                    seq_data = data[seq_start_rel:seq_end_rel]
                    
                    if len(seq_data) == 0:
                        cursor = seq_end_rel
                        continue
                    
                    cursor = seq_end_rel
                    
                else:
                    header_start_rel = data.find(b'@', cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                    if header_start_rel == -1:
                        break
                    
                    header_end_rel = data.find(b'\n', header_start_rel, min(len(data), header_start_rel + 10000))
                    if header_end_rel == -1:
                        break
                    
                    header = data[header_start_rel:header_end_rel].decode('utf-8', errors='ignore')
                    
                    seq_start_rel = header_end_rel + 1
                    
                    next_at_rel = data.find(b'@', seq_start_rel + 1, min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION))
                    
                    if next_at_rel != -1 and next_at_rel > seq_start_rel:
                        if data[next_at_rel-1:next_at_rel] == b'\n':
                            seq_end_rel = next_at_rel - 1
                        else:
                            seq_end_rel = next_at_rel
                    else:
                        seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                        if seq_end_rel > seq_start_rel and data[seq_start_rel:seq_end_rel][-1:] == b'\n':
                            seq_end_rel -= 1
                    
                    seq_data = data[seq_start_rel:seq_end_rel]
                    cursor = seq_end_rel + 1

                # Check for metadata block changes (common for both safe and non-safe mode)
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_inverse_table = inverse_tables[current_metadata_idx + 1]
                        current_metadata_idx += 1

                # Convert binary to bases
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
                
                # Write output
                output_header = f"{header}\n{bases}\n+\n" if not second_head_flag else f"{header}\n{bases}\n+{header[1:]}\n"
                output_buffer.write(output_header.encode('ascii'))
                output_buffer.write(quality_string.encode('ascii'))
                output_buffer.write(b'\n')
                sequence_count += 1
                sequences_in_chunk += 1
            
            return (chunk_id, output_buffer.getvalue(), sequences_in_chunk)

        
        elif mode == 0:
            # Mode 0: Headers compressed, bases and quality as plain text 
            # @header\nBASES\n+\nQUALITY\n
            sequences_in_chunk = 0
            cursor = 0
            chunk_size = abs_end - abs_start
            
            while cursor < chunk_size:
                if data[cursor:cursor+1] != b'@':
                    # Try to find next @ 
                    next_at = data.find(b'@', cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                    if next_at == -1:
                        break
                    cursor = next_at
                
                line1_end = data.find(b'\n', cursor) # Header
                if line1_end == -1 or line1_end >= len(data):
                    break
                
                compressed_header = data[cursor+1:line1_end].decode('utf-8', errors='replace').strip()
                pair_number = 0
                
                line2_start = line1_end + 1
                seq_cursor = line2_start
                bases_len = 0
                line3_start = -1
                bases_parts = []
                
                while seq_cursor < len(data):
                    line_end = data.find(b'\n', seq_cursor)
                    if line_end == -1:
                        break
                    
                    next_line_start = line_end + 1
                    if next_line_start < len(data) and data[next_line_start:next_line_start+1] == b'+':
                        bases_parts.append(data[seq_cursor:line_end])
                        bases_len = sum(len(p) for p in bases_parts)
                        line3_start = next_line_start
                        break
                    
                    bases_parts.append(data[seq_cursor:line_end])
                    seq_cursor = line_end + 1
                
                if line3_start == -1:
                    break
                
                bases_line = b''.join(bases_parts).decode('utf-8', errors='replace').strip()
                line3_end = data.find(b'\n', line3_start) # Line 3
                if line3_end == -1:
                    break
                
                plus_line = data[line3_start:line3_end].decode('utf-8', errors='replace').strip()
                
                if not plus_line.startswith('+'):
                    cursor = line3_end + 1
                    continue
                
                line4_start = line3_end + 1
                line4_end = line4_start + bases_len
                
                if line4_end > len(data):
                    break
                
                quality_line = data[line4_start:line4_end].decode('utf-8', errors='replace')
                
                unique_id = compressed_header
                
                # Check for metadata block changes
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_metadata_idx += 1
                
                # Reconstruct full header
                if current_metadata and current_metadata.structure_template:
                    full_header = reconstruct_header_from_structure(
                        current_metadata.structure_template,
                        unique_id,
                        current_metadata.sequencer_type,
                        pair_number
                    )
                else:
                    full_header = f"@{unique_id}"
                    if pair_number > 0:
                        full_header = f"{full_header}/{pair_number}"
                
                # Add length flag if needed
                if length_flag:
                    full_header = f"{full_header} length={len(bases_line)}"
                
                # Write reconstructed record
                output_buffer.write(full_header.encode('utf-8'))
                output_buffer.write(b'\n')
                output_buffer.write(bases_line.encode('utf-8'))
                output_buffer.write(b'\n')
                
                if second_head_flag:
                    output_buffer.write(b'+')
                    output_buffer.write(full_header[1:].encode('utf-8'))  # Remove '@'
                    output_buffer.write(b'\n')
                else:
                    output_buffer.write(b'+\n')
                
                output_buffer.write(quality_line.encode('utf-8'))
                output_buffer.write(b'\n')
                
                sequence_count += 1
                sequences_in_chunk += 1
                
                # Move cursor past this record (skip trailing newline if present)
                cursor = line4_end
                if cursor < len(data) and data[cursor:cursor+1] == b'\n':
                    cursor += 1
            
            return (chunk_id, output_buffer.getvalue(), sequences_in_chunk)
        

        elif mode == 3:
            sequences_in_chunk = 0
            cursor = 0
            chunk_size = abs_end - abs_start
            
            MAX_CHUNK_EXTENSION = 10 * 1024 * 1024  
            
            while cursor < chunk_size:
                marker_pos = data.find(b'\xff', cursor)
                
                if marker_pos == -1:
                    extension_needed = min(MAX_CHUNK_EXTENSION, len(data) - cursor)
                    if cursor + extension_needed >= len(data):
                        break
                    marker_pos = data.find(b'\xff', cursor, cursor + extension_needed)
                    if marker_pos == -1:
                        break
                
                seq_start_rel = marker_pos + 1
                next_marker = data.find(b'\xff', seq_start_rel)
                
                if next_marker != -1:
                    seq_end_rel = next_marker
                else:
                    lookahead_end = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                    seq_end_rel = data.find(b'\xff', seq_start_rel, lookahead_end)
                    if seq_end_rel == -1:
                        seq_end_rel = lookahead_end
                
                seq_data = data[seq_start_rel:seq_end_rel]
                
                if len(seq_data) > 0 and seq_data[-1:] == b'\n':
                    seq_data = seq_data[:-1]
                
                if len(seq_data) == 0:
                    cursor = marker_pos + 1
                    continue
                
                # Check for metadata block changes
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_inverse_table = inverse_tables[current_metadata_idx + 1]
                        current_metadata_idx += 1

                # Convert binary to bases
                seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                bases_array = reverse_map[seq_array]
                bases = bases_array.tobytes().decode('ascii')
                
                if header_index is not None and sequence_count < len(header_index):
                    header = get_header_by_index(headers_file_path, header_index, sequence_count)
                    if header and not header.startswith('@'):
                        header = '@' + header
                else:
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

                # Write output
                output_buffer.write(header.encode('utf-8'))
                output_buffer.write(b'\n')
                output_buffer.write(bases.encode('ascii'))
                
                plus_line = b"\n+"
                if second_head_flag and header:
                    plus_line += header[1:].encode('utf-8')
                plus_line += b"\n"
                output_buffer.write(plus_line)
                
                output_buffer.write(quality_string.encode('ascii'))
                output_buffer.write(b'\n')
                
                sequence_count += 1 
                sequences_in_chunk += 1
                local_count += 1
                cursor = seq_end_rel 
            
            return (chunk_id, output_buffer.getvalue(), local_count)

        elif mode == 2:
            # Safe mode, format is: \n@header(255)\nbases\n
            sequences_in_chunk = 0
            cursor = 0
            chunk_size = abs_end - abs_start
            
            while cursor < chunk_size:
                if safe_mode_flag:
                    # Find next @ with validation
                    if cursor == 0 and data[0:1] == b'@':
                        # First sequence in chunk starts at position 0
                        header_start_rel = 0
                    else:
                        # Always look for \n@ for subsequent sequences
                        header_start_rel = data.find(b'\n@', cursor)
                        if header_start_rel != -1:
                            header_start_rel += 1 # Skip the \n, point to the @
                    
                    if header_start_rel < 0 or header_start_rel >= chunk_size:
                        break
                    
                    candidate = header_start_rel # Then validate by finding \xff after this
                    
                    while candidate >= 0:
                        next_xff_rel = data.find(b'\xff', candidate, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                        next_at_rel = data.find(b'\n@', candidate + 1, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                        
                        # If another \n@ comes before \xff, current candidate is invalid
                        if next_at_rel != -1 and (next_xff_rel == -1 or next_at_rel < next_xff_rel):
                            candidate = next_at_rel + 1
                            if candidate >= chunk_size + MAX_CHUNK_EXTENSION:
                                candidate = -1
                                break
                            continue
                        
                        # Found valid \xff
                        if next_xff_rel != -1:
                            xff_pos = next_xff_rel
                            break
                        
                        candidate = -1
                        break
                    
                    if candidate == -1:
                        break
                    
                    header_start = candidate
                    
                    # Extract header (between @ and \xff) w/ mode 2 reconstructing the header based off structure
                    header_content = data[header_start + 1:xff_pos].decode('utf-8', errors='replace').strip()
                    
                    # Bases start after \xff\n or \xff
                    if xff_pos + 1 < len(data) and data[xff_pos+1:xff_pos+2] == b'\n':
                        seq_start_rel = xff_pos + 2
                    else:
                        seq_start_rel = xff_pos + 1
                    
                    # Find sequence end (next valid @ header)
                    candidate_end_rel = seq_start_rel
                    seq_end_rel = -1
                    
                    while True:
                        next_at_rel = data.find(b'\n@', candidate_end_rel, min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION))
                        
                        if next_at_rel != -1:
                            # Validate our \n@
                            next_xff_rel = data.find(b'\xff', next_at_rel + 1, min(len(data), next_at_rel + 1000))
                            following_at_rel = data.find(b'\n@', next_at_rel + 2, min(len(data), next_at_rel + 1000))
                            
                            if following_at_rel != -1 and (next_xff_rel == -1 or following_at_rel < next_xff_rel):
                                candidate_end_rel = following_at_rel + 1
                                if candidate_end_rel > seq_start_rel + MAX_CHUNK_EXTENSION:
                                    break
                                continue
                            
                            # Valid header if sequence ends at the \n
                            seq_end_rel = next_at_rel
                            break
                        else:
                            seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                            while seq_end_rel > seq_start_rel and data[seq_end_rel-1:seq_end_rel] in (b'\n', b'\r', b' '):
                                seq_end_rel -= 1
                            break
                    
                    if seq_end_rel == -1 or seq_end_rel <= seq_start_rel:
                        cursor = len(data)
                        continue
                    
                    seq_data = data[seq_start_rel:seq_end_rel]
                    
                    if len(seq_data) == 0:
                        cursor = seq_end_rel
                        continue
                    
                    cursor = seq_end_rel
                    
                else:
                    # Search in full file
                    start_idx_rel = data.find(b'@', cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION))
                    if start_idx_rel == -1:
                        break
                    
                    header_end_rel = data.find(b'\n', start_idx_rel, min(len(data), start_idx_rel + 10000))
                    if header_end_rel == -1:
                        break
                    
                    header_content = data[start_idx_rel+1:header_end_rel].decode('utf-8', errors='ignore').strip()
                    
                    seq_start_rel = header_end_rel + 1
                    
                    next_header_rel = data.find(b'\n@', seq_start_rel, min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION))
                    
                    if next_header_rel != -1:
                        seq_end_rel = next_header_rel
                    else:
                        seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                        while seq_end_rel > seq_start_rel and data[seq_end_rel-1:seq_end_rel] in (b'\n', b'\r', b' ', b'\t'):
                            seq_end_rel -= 1
                    
                    seq_data = data[seq_start_rel:seq_end_rel]
                    cursor = seq_end_rel + 1
                
                unique_id = None
                pair_number = 0
                
                if '/' in header_content:
                    parts = header_content.rsplit('/', 1)
                    if len(parts) == 2 and parts[1].isdigit() and len(parts[1]) == 1 and int(parts[1]) in [1, 2]:
                        unique_id = parts[0]
                        pair_number = int(parts[1])
                    else:
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
                        header = f"@{sra_accession}.{sequence_count}"
                    else:
                        header = f"@seq{sequence_count}"
                
                if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                    next_metadata = metadata_blocks[current_metadata_idx + 1]
                    if sequence_count >= next_metadata.start_index:
                        current_metadata = next_metadata
                        current_inverse_table = inverse_tables[current_metadata_idx + 1]
                        current_metadata_idx += 1
                
                seq_array = np.frombuffer(seq_data, dtype=np.uint8)
                bases_array = reverse_map[seq_array]
                bases = bases_array.tobytes().decode('ascii')
                
                if current_inverse_table is not None:
                    quality_scores = reverse_scaling_to_quality(
                        seq_array, subtract_table,
                        current_inverse_table, max_phred=phred_alphabet_max
                    )
                    quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
                
                if length_flag:
                    header = f"{header} length={len(bases)}"
                
                output_buffer.write(header.encode('utf-8'))
                output_buffer.write(b'\n')
                output_buffer.write(bases.encode('ascii'))
                output_buffer.write(b'\n+')
                if second_head_flag:
                    output_buffer.write(header[1:].encode('utf-8'))
                output_buffer.write(b'\n')
                output_buffer.write(quality_string.encode('ascii'))
                output_buffer.write(b'\n')
                
                sequence_count += 1
                sequences_in_chunk += 1
            
            return (chunk_id, output_buffer.getvalue(), sequences_in_chunk)

    except Exception as e:
        logger.warning(f"ERROR in worker processing chunk {chunk_id}: {e}")
        traceback.print_exc()
        raise


def reconstruct_fastq(input_path: str, output_path: str, 
                     gray_N: int = 0, gray_A: int = 3, gray_G: int = 66,
                     gray_C: int = 129, gray_T: int = 192, phred_alphabet_max: int = None,
                     phred_offset: int = 33, chunk_size_mb: int = 8, num_workers: int = 4,
                     mode: int = 2, mode3_headers_file: str = None, verbose: bool=False):
    """
    Reconstruct FASTQ file from FASTR.
    """
    logger.info(f"Reading FASTR: {input_path}")
    logger.info(f"Reconstruction mode: {mode}")
    
    headers_mmap_info = None
    header_index_mmap = None
    mmap_cleanup_path = None
    
    if mode == 3 and mode3_headers_file:
        if not os.path.exists(mode3_headers_file):
            raise FileNotFoundError(f"Headers file not found: {mode3_headers_file}")
        
        header_index_mmap, mmap_cleanup_path = build_header_index(mode3_headers_file)
        headers_mmap_info = (mmap_cleanup_path, header_index_mmap.shape, header_index_mmap.dtype)
        logger.info(f"Header index created with {len(header_index_mmap):,} entries")
        
    elif mode == 3 and not mode3_headers_file:
        logger.warning("WARNING: Mode 3 requires --headers_file argument")
        logger.info("Proceeding without headers - will use fallback header generation")

    file_size = os.path.getsize(input_path)
    logger.info(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    
    logger.info(f"Using {num_workers} worker processes for parallel reconstruction")
    
    HEADER_READ_SIZE = 10 * 1024 * 1024  
    with open(input_path, 'rb') as f:
        header_data = f.read(HEADER_READ_SIZE)
    
    metadata_blocks, data_start_byte, sra_accession, phred_from_metadata, detected_mode, length_flag, second_head_flag, safe_mode_flag = parse_metadata_header(header_data, mode)

    if detected_mode is not None:
        mode = detected_mode
        logger.info(f"Detected mode from metadata: {mode}")
    
    logger.info(f"Reconstruction mode: {mode}")
    
    if phred_alphabet_max is None:
        if phred_from_metadata is not None:
            phred_alphabet_max = phred_from_metadata
            logger.info(f"Using PHRED alphabet from metadata: {phred_alphabet_max}")
        else:
            phred_alphabet_max = 41
            logger.info(f"No PHRED alphabet found in metadata, using default: {phred_alphabet_max}")
    else:
        logger.info(f"Using user-specified PHRED alphabet (overriding metadata): {phred_alphabet_max}")
    
    if mode in [0, 2] and not metadata_blocks:
        logger.warning("WARNING: No metadata found for header reconstruction")
    
    if metadata_blocks:
        logger.info(f"Found {len(metadata_blocks)} metadata block(s)")
        for i, mb in enumerate(metadata_blocks):
            logger.info(f"  Block {i+1}: {mb.sequencer_type}, equation: {mb.scaling_equation}")
            if mb.structure_template:
                logger.info(f"           structure: {mb.structure_template}")
    
    if sra_accession:
        logger.info(f"SRA Accession: {sra_accession}")
    
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
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    logger.info(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    def chunk_generator():
        chunk_id = 0
        start_seq_idx = 0
        
        if mode == 0:
            logger.info("Mode 0: Building record index...")
            record_positions = [data_start_byte]
            
            with open(input_path, 'rb') as f:
                f.seek(data_start_byte)
                SCAN_BLOCK_SIZE = 50 * 1024 * 1024
                buffer = b''
                file_pos = data_start_byte
                
                while True:
                    chunk = f.read(SCAN_BLOCK_SIZE)
                    if not chunk:
                        break
                    
                    buffer += chunk
                    cursor = 0
                    # We run off the assumption each fastq record is exactly 4 lines, so we'll count...
                    while cursor < len(buffer):
                        
                        if buffer[cursor:cursor+1] != b'@':
                            cursor += 1
                            continue
                        line1_end = buffer.find(b'\n', cursor)
                        if line1_end == -1:
                            break
                        
                        line2_start = line1_end + 1
                        line2_end = buffer.find(b'\n', line2_start)
                        if line2_end == -1:
                            break
                        
                        bases_len = line2_end - line2_start
                        line3_start = line2_end + 1
                        line3_end = buffer.find(b'\n', line3_start)
                        if line3_end == -1:
                            break
                        
                        if buffer[line3_start:line3_start+1] != b'+':
                            cursor += 1
                            continue
                        
                        line4_start = line3_end + 1
                        line4_end = line4_start + bases_len
                        
                        if line4_end > len(buffer):
                            break
                        
                        cursor = line4_end
                        if cursor < len(buffer) and buffer[cursor:cursor+1] == b'\n':
                            cursor += 1
                        if cursor < len(buffer):
                            record_positions.append(file_pos + cursor)
                    
                    buffer = buffer[cursor:]
                    file_pos += cursor
            
            logger.info(f"Mode 0: Found {len(record_positions)} record boundaries")
            
            target_bytes_per_chunk = chunk_size_bytes
            current_chunk_start = 0
            current_chunk_byte_count = 0
            
            for i in range(len(record_positions)):
                if i == 0:
                    continue
                
                record_start = record_positions[i-1]
                record_end = record_positions[i] if i < len(record_positions) - 1 else file_size
                record_bytes = record_end - record_start
                
                if (current_chunk_byte_count + record_bytes > target_bytes_per_chunk and 
                    current_chunk_byte_count > 0): 
                    
                    abs_start = record_positions[current_chunk_start]
                    abs_end = record_positions[i-1]
                    seqs_in_chunk = (i-1) - current_chunk_start
                    
                    if verbose:
                        logger.info(f"Yielding chunk {chunk_id}: {abs_start}-{abs_end} " +
                                f"({abs_end-abs_start:,} bytes, {seqs_in_chunk} seqs)")
                    
                    yield (chunk_id, abs_start, abs_end, start_seq_idx) # If adding this record would exceed chunk size, yield current chunk
                    
                    start_seq_idx += seqs_in_chunk
                    chunk_id += 1
                    current_chunk_start = i-1
                    current_chunk_byte_count = record_bytes
                else:
                    current_chunk_byte_count += record_bytes
            
            if current_chunk_start < len(record_positions):
                abs_start = record_positions[current_chunk_start]
                abs_end = file_size
                seqs_in_chunk = len(record_positions) - current_chunk_start
                
                if verbose:
                    logger.info(f"Yielding final chunk {chunk_id}: {abs_start}-{abs_end} " +
                            f"({abs_end-abs_start:,} bytes, {seqs_in_chunk} seqs)")
                
                yield (chunk_id, abs_start, abs_end, start_seq_idx) # Yield the final chunk
            
            return
        
        # For modes 1, 2, 3: Stream through file finding chunk boundaries
        buffer_start = data_start_byte
        buffer_len = file_size - data_start_byte
        pos = 0
        
        with open(input_path, 'rb') as f:
            while pos < buffer_len:
                chunk_end = min(pos + chunk_size_bytes, buffer_len)
                
                # Read search region for boundary detection
                search_start = max(pos, chunk_end - 1000)
                search_size = chunk_end - search_start + 1000  # Extra buffer for safety
                
                f.seek(buffer_start + search_start)
                search_region = f.read(min(search_size, buffer_len - search_start))
                
                # Find boundary based on mode
                if mode == 3:
                    # Mode 3 splits on 255 markers
                    last_marker = search_region.rfind(b'\n\xff', 0, chunk_end - search_start)
                    
                    if last_marker != -1:
                        chunk_end = search_start + last_marker + 1
                
                elif mode == 2 and safe_mode_flag:
                    # Safe mode splits on validated @header\xff pattern
                    candidate_pos = min(len(search_region) - 1, chunk_end - search_start)
                    found_valid = False
                    
                    while candidate_pos > 0:
                        at_pos = search_region.rfind(b'@', 0, candidate_pos)
                        if at_pos == -1:
                            break
                        
                        # Check if this @ has \xff after it (within reasonable distance)
                        xff_pos = search_region.find(b'\xff', at_pos, min(at_pos + 200, len(search_region)))
                        intervening_at = search_region.find(b'\n@', at_pos + 1, min(at_pos + 200, len(search_region)))
                        
                        # Valid only if we found \xff AND no \n@ comes before it
                        if xff_pos != -1 and (intervening_at == -1 or xff_pos < intervening_at):
                            chunk_end = search_start + at_pos
                            found_valid = True
                            break
                        
                        candidate_pos = at_pos
                    
                    if not found_valid:
                        chunk_end = min(pos + chunk_size_bytes, buffer_len)
                
                else:
                    # Mode 1 and 2 (unsafe) split on \n@ headers
                    last_header = search_region.rfind(b'\n@', 0, chunk_end - search_start)
                    
                    if last_header != -1:
                        chunk_end = search_start + last_header + 1
                
                if chunk_end <= pos:
                    chunk_end = min(pos + chunk_size_bytes, buffer_len)
                
                abs_start = buffer_start + pos
                abs_end = buffer_start + chunk_end
                
                if abs_end > abs_start:
                    # Estimate sequences for this chunk
                    f.seek(abs_start)
                    sample = f.read(min(abs_end - abs_start, 100000))  # Sample first 100KB
                    
                    if mode == 3:
                        estimated_seqs = sample.count(b'\xff')
                        if len(sample) < abs_end - abs_start:
                            # Extrapolate
                            estimated_seqs = int(estimated_seqs * (abs_end - abs_start) / len(sample))
                    else:
                        estimated_seqs = sample.count(b'\n@')
                        if sample.startswith(b'@'):
                            estimated_seqs += 1
                        if len(sample) < abs_end - abs_start:
                            estimated_seqs = int(estimated_seqs * (abs_end - abs_start) / len(sample))
                    
                    if verbose:
                        logger.info(f"Yielding chunk {chunk_id}: {abs_start}-{abs_end} ({abs_end-abs_start} bytes, ~{estimated_seqs} seqs)")
                    
                    yield (chunk_id, abs_start, abs_end, start_seq_idx)
                    
                    start_seq_idx += estimated_seqs
                    chunk_id += 1
                
                pos = chunk_end
                if pos >= buffer_len:
                    break

    logger.info("Reconstructing FASTQ...")

    chunks_list = []
    current_seq_idx = 0

    for chunk_data in chunk_generator():
        chunk_id, abs_start, abs_end, _ = chunk_data 
        
        with open(input_path, 'rb') as f:
            f.seek(abs_start)
            chunk_size = abs_end - abs_start
            chunk_bytes = f.read(chunk_size)
            actual_count = chunk_bytes.count(b'\xff') 
        
        chunks_list.append((chunk_id, abs_start, abs_end, current_seq_idx))
        current_seq_idx += actual_count

    total_sequences = 0

    try:
        with open(output_path, 'wb', buffering=chunk_size_bytes) as outfile:
            with Pool(processes=num_workers) as pool:
                worker_func = partial( 
                    process_chunk_worker_reconstruction,
                    mmap_path=input_path,        
                    reverse_map=reverse_map,       
                    subtract_table=subtract_table,  
                    metadata_blocks=metadata_blocks,
                    inverse_tables=inverse_tables,
                    phred_alphabet_max=phred_alphabet_max,
                    phred_offset=phred_offset,
                    sra_accession=sra_accession,
                    mode=mode,
                    headers_file_path=mode3_headers_file,
                    headers_mmap_info=headers_mmap_info, 
                    data_start_byte=data_start_byte,
                    length_flag=length_flag,
                    second_head_flag=second_head_flag,
                    safe_mode_flag=safe_mode_flag
                )
                
                for chunk_id, processed_bytes, count in pool.imap(worker_func, chunks_list, chunksize=1):
                    outfile.write(processed_bytes)
                    total_sequences += count
    
    finally:
        if header_index_mmap is not None:
            del header_index_mmap  # Close mmap
        if mmap_cleanup_path and os.path.exists(mmap_cleanup_path):
            try:
                os.unlink(mmap_cleanup_path)
                logger.info("Cleaned up temporary index file")
            except:
                pass
    
    logger.info(f"\nTotal sequences reconstructed: {total_sequences:,}")
    logger.info(f"Output saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct FASTQ files from FASTR.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Positional Arguments
    parser.add_argument("input_path",
                        metavar="FILE",
                        help="Path to FASTR compressed file")
    parser.add_argument("output_path", 
                        metavar="FILE",
                        help="Output FASTQ file path")

    # Mode Group
    mode_group = parser.add_argument_group("RECONSTRUCTION MODE")
    mode_group.add_argument("--mode", type=int, metavar="INT", default=2,
                            choices=[0, 1, 2, 3],
                            help="Reconstruction mode [2]\n"
                                 "0: Headers only (no base conversion)\n"
                                 "1: Bases only (keep original headers)\n"
                                 "2: Full reconstruction (headers + bases)\n"
                                 "3: No repeating headers (requires --headers_file)")
    mode_group.add_argument("--headers_file", type=str, metavar="FILE", default=None,
                            help="Path to headers file for mode 3 reconstruction [null]")

    # Quality Parameters Group
    quality_group = parser.add_argument_group("QUALITY RECONSTRUCTION")
    quality_group.add_argument("--phred_offset", type=int, metavar="INT", default=33,
                               help="Phred quality offset for output [33]")
    quality_group.add_argument("--phred_alphabet", type=str, metavar="STR", default=None,
                               choices=['phred42', 'phred63', 'phred94'],
                               help="Override phred alphabet from metadata (phred42/phred63/phred94) [auto]")

    # Base Mapping Group
    gray_group = parser.add_argument_group("GRAYSCALE DECODING")
    gray_group.add_argument("--gray_N", type=int, metavar="INT", default=0,
                            help="Grayscale value for N [0]")
    gray_group.add_argument("--gray_A", type=int, metavar="INT", default=3,
                            help="Grayscale value for A [3]")
    gray_group.add_argument("--gray_G", type=int, metavar="INT", default=66,
                            help="Grayscale value for G [66]")
    gray_group.add_argument("--gray_C", type=int, metavar="INT", default=129,
                            help="Grayscale value for C [129]")
    gray_group.add_argument("--gray_T", type=int, metavar="INT", default=192,
                            help="Grayscale value for T [192]")

    # Performance Group
    perf_group = parser.add_argument_group("PERFORMANCE & PARALLELIZATION")
    perf_group.add_argument("--chunk_size_mb", type=int, metavar="INT", default=8,
                            help="Chunk size in MB for parallel processing [8]")
    perf_group.add_argument("--num_workers", type=int, metavar="INT", default=4,
                            help="Number of parallel workers [4]")
    perf_group.add_argument("--verbose", type=int, metavar="INT", default=0,
                            choices=[0, 1],
                            help="Enable verbose logging (0/1) [0]")
    perf_group.add_argument("--profile", type=int, metavar="INT", default=0,
                            choices=[0, 1],
                            help="Enable cProfile profiling (0/1) [0]")

    args = parser.parse_args()

    if args.verbose == 1:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    phred_alphabet_max = None
    if args.phred_alphabet:
        if args.phred_alphabet == "phred42":
            phred_alphabet_max = 41
        elif args.phred_alphabet == "phred63":
            phred_alphabet_max = 62
        elif args.phred_alphabet == "phred94":
            phred_alphabet_max = 93
    
    start_time = time.perf_counter()

    if args.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        logger.info("Profiling enabled...")
    
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
        mode3_headers_file=args.headers_file,
        verbose=(args.verbose == 1)
    )
    if args.profile == 1:
        profiler.disable()
        s = StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
        ps.print_stats(20)
        print("\n" + "="*80)
        print("Profiling Results:")
        print("="*80)
        print(s.getvalue())
    
    end_time = time.perf_counter()
    logger.info(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()