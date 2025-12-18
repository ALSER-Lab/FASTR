import os
import numpy as np
import time
import argparse
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

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
    
    first_seq_marker = data.find(b'<') # Find where actual sequence data starts (which is marked by a '<')
    if first_seq_marker == -1:
        return metadata_blocks, 0, sra_accession
    
    header_section = data[:first_seq_marker].decode('utf-8', errors='ignore') # Extract ONLY the header section before the first '<'
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
    
    if not metadata_blocks:
        print("WARNING: No metadata blocks found, using defaults") 
        # Create a default block if none found
        metadata_blocks.append(MetadataBlock(
            common_metadata={'prefix': sra_accession.split('.')[0][:3] if sra_accession else '', 
                           'accession': sra_accession.split('.')[0][3:] if sra_accession else ''} if sra_accession else {},
            sequencer_type='srr' if sra_accession else 'none',
            scaling_equation='x',
            start_index=0,
            end_index=-1
        ))
    
    return metadata_blocks, first_seq_marker, sra_accession

def parse_custom_formula(formula: str, quality_scores: np.ndarray) -> np.ndarray:
    # Just copy-pasted from sequence_converter.py
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
    for i in range(1, 64):
        reverse_map[i] = ord('A')
    
    # T: 65-127
    for i in range(65, 128):
        reverse_map[i] = ord('T')

    # C: 129-191
    for i in range(129, 192):
        reverse_map[i] = ord('C')
    
    # G: 193-255
    for i in range(193, 256):
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
    
    # Clip to valid range BEFORE rounding
    scaled_for_q = np.clip(scaled_for_q, 0, max_range)
    
    # Round to integers
    scaled_int_for_q = np.round(scaled_for_q).astype(int)
    
    # Build inverse lookup table for scaled values 0..max_range
    inverse_table = np.full(max_range + 1, -1, dtype=int)
    
    # For each quality score, store it at its corresponding scaled position
    for q, s in zip(q_possible, scaled_int_for_q):
        if 0 <= s <= max_range:
            inverse_table[s] = int(q)
    
    # Fill gaps in inverse table using interpolation
    valid_indices = np.where(inverse_table != -1)[0]
    
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


def reconstruct_fastq(input_path: str, output_path: str, 
                     gray_N: int = 1, gray_A: int = 63, gray_T: int = 127,
                     gray_C: int = 191, gray_G: int = 255, phred_alphabet_max: int = 41,
                     phred_offset: int = 33):
    """
    Reconstruct FASTQ file from FASTR (mode 3).
    Assumes binary writing was used in sequence_converter (--binary_write 1).
    """
    print(f"Reading compressed binary file: {input_path}")
    
    # Read entire file as binary (this can vary, and I will add auto-detection to sequence_converter.py to determine if its written in binary/ints),
    # But binary is what we will be testing with (as it is more important, and more likely to be used)
    with open(input_path, 'rb') as f:
        data = f.read()
    file_size = len(data)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**2):.2f} MB)")
    metadata_blocks, data_start_byte, sra_accession = parse_metadata_header(data) # Parse metadata from text header
    
    if not metadata_blocks:
        print("ERROR: No metadata found?")
        return
    
    print(f"Found {len(metadata_blocks)} metadata block(s)")
    if sra_accession:
        print(f"SRA Accession: {sra_accession}")
    
    for i, mb in enumerate(metadata_blocks):
        print(f"  Block {i+1}: {mb.sequencer_type}, equation: {mb.scaling_equation}")
    
    # Create reverse base mapping
    reverse_map = reverse_base_map(gray_N, gray_A, gray_T, gray_C, gray_G)
    
    # Pre-compute inverse quality tables for each metadata block
    inverse_tables = []
    for mb in metadata_blocks:
        formula_func = build_formula_func(mb.scaling_equation)
        inverse_table = build_inverse_quality_table(formula_func, phred_alphabet_max)
        inverse_tables.append(inverse_table)
    
    # Base ranges for quality decoding
    base_ranges = {
        'A': (1, 63),
        'T': (65, 127),
        'C': (129, 191),
        'G': (193, 255),
        'N': (1, 1)
    }
    
    # Process binary sequence data
    print(f"Reconstructing FASTQ from byte {data_start_byte}...")
    
    sequence_count = 0
    current_metadata_idx = 0
    current_metadata = metadata_blocks[0]
    current_inverse_table = inverse_tables[0]
    
    # Calculate 10% milestones based on file size
    data_size = file_size - data_start_byte
    milestone_interval = data_size // 10
    next_milestone = milestone_interval
    last_reported_percent = 0
    output_buffer = []
    buffer_size = 10000
    
    # Start reading from after metadata
    pos = data_start_byte
    with open(output_path, 'w') as outfile:
        while pos < len(data):
            # Check for 10% progress milestones
            bytes_processed = pos - data_start_byte
            if milestone_interval > 0 and bytes_processed >= next_milestone:
                percent_complete = int((bytes_processed / data_size) * 100)
                if percent_complete > last_reported_percent:
                    print(f"Progress: {percent_complete}% complete ({sequence_count:,} sequences)")
                    last_reported_percent = percent_complete
                next_milestone += milestone_interval
            
            # Find sequence start marker, which is denoted by a '<'
            if data[pos:pos+1] != b'<':
                pos += 1
                continue
            
            pos += 1  # Skip '<'
            
            # Find newline to get sequence length
            seq_end = data.find(b'\n', pos)
            if seq_end == -1:
                # Check if there's still data left (last sequence might not have trailing newline)
                if pos < len(data):
                    seq_end = len(data)
                else:
                    break
            
            # Extract binary sequence data
            seq_data = data[pos:seq_end]
            
            # Check if we need to switch metadata blocks (different flowcells will be given different metadata blocks, though this is likely to be subject to change)
            if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
                next_metadata = metadata_blocks[current_metadata_idx + 1]
                if sequence_count >= next_metadata.start_index:
                    current_metadata = next_metadata
                    current_inverse_table = inverse_tables[current_metadata_idx + 1]
                    current_metadata_idx += 1
                    print(f"Switched to metadata block {current_metadata_idx + 1}")
            
            # Convert binary values to bases
            seq_array = np.frombuffer(seq_data, dtype=np.uint8)
            bases_array = reverse_map[seq_array]
            bases = bases_array.tobytes().decode('ascii')
            
            # Try to find header before this sequence
            # Look backwards for '@' followed by text and newline
            header_search_start = max(0, pos - 500)  # Search up to 500 bytes back
            header_section = data[header_search_start:pos-1]
            
            unique_id = None
            pair_number = 0
            
            # Find last '@' before '<'
            last_at = header_section.rfind(b'@')
            if last_at != -1:
                # Extract potential header
                header_line_start = header_search_start + last_at
                header_line_end = data.find(b'\n', header_line_start)
                if header_line_end != -1 and header_line_end < pos:
                    header_content = data[header_line_start+1:header_line_end].decode('utf-8', errors='ignore').strip()
                    
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
            
            # Reconstruct header
            if unique_id:
                header = reconstruct_header(unique_id, current_metadata.common_metadata,
                                            current_metadata.sequencer_type, pair_number, sra_accession)
            else:
                # Fallback header with SRA accession if available
                if sra_accession:
                    header = f"@{sra_accession}.{sequence_count}"
                else:
                    header = f"@seq_{sequence_count}"
            
            # Apply reverse scaling using pre-computed inverse table
            quality_scores = reverse_scaling_to_quality(
                seq_array,
                base_ranges,
                reverse_map,
                current_inverse_table,
                max_phred=phred_alphabet_max
            )
            quality_string = quality_to_ascii(quality_scores, phred_offset).decode('ascii')
            
            # Add to output buffer
            output_buffer.append(f"{header}\n{bases}\n+\n{quality_string}\n")
            
            # Write buffer when it reaches threshold
            if len(output_buffer) >= buffer_size:
                outfile.write(''.join(output_buffer))
                output_buffer.clear()
            
            sequence_count += 1
            # Move to next sequence
            pos = seq_end + 1
        
        # Write any remaining buffered output
        if output_buffer:
            outfile.write(''.join(output_buffer))
    
    print(f"Progress: 100% complete")
    print(f"\nTotal sequences reconstructed: {sequence_count:,}")
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
        phred_offset=args.phred_offset
    )
    
    end_time = time.perf_counter()
    print(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()