import numpy as np
from typing import List, Tuple, Optional, Dict
from data_structures import FASTQRecord
from header_compression import compress_header


def parse_fastq_records_from_buffer(buffer: bytes, start_index: int, base_map: np.ndarray,
                                    phred_map: Optional[np.ndarray], compress_headers: bool,
                                    sequencer_type: str, paired_end: bool, keep_bases: bool,
                                    binary_bases: bool, keep_quality: bool, binary_quality: bool) -> Tuple[List[FASTQRecord], bytes, Optional[Dict], int]:
    """
    Parse FASTQ records from a buffer using vectorized operations.
    Returns (records, leftover_buffer, current_flowcell_metadata, num_records)
    """
    
    records = []
    current_flowcell_metadata = None
    lines = buffer.split(b'\n')
    
    # Check if we have complete 4-line records
    if lines and lines[-1] == b'':
        lines = lines[:-1]
    
    num_complete_records = len(lines) // 4
    
    if num_complete_records == 0:
        return records, buffer, current_flowcell_metadata, 0
    
    leftover = b'\n'.join(lines[num_complete_records * 4:])
    
    # Pre-create binary map ONCE if needed
    binary_map = None
    if binary_bases:
        binary_map = np.zeros(128, dtype=np.uint8)
        binary_map[ord('A')] = 0
        binary_map[ord('T')] = 1
        binary_map[ord('C')] = 2
        binary_map[ord('G')] = 3
        binary_map[ord('N')] = 4
    
    # Process in batches to avoid numpy split overhead
    BATCH_SIZE = 10000  # Process 10k records at a time
    
    for batch_start in range(0, num_complete_records, BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, num_complete_records)
        batch_size = batch_end - batch_start
        
        line_start = batch_start * 4
        line_end = batch_end * 4
        
        batch_lines = lines[line_start:line_end]
        
        # Vectorized operations on this batch
        header_strs = [batch_lines[i].decode('utf-8', errors='ignore') for i in range(0, len(batch_lines), 4)]
        seq_data_list = [batch_lines[i].replace(b'\r', b'') for i in range(1, len(batch_lines), 4)]
        qual_data_list = [batch_lines[i].replace(b'\r', b'') for i in range(3, len(batch_lines), 4)]
        
        # Concatenate and convert ONCE for entire batch (major speedup)
        seq_lengths = [len(s) for s in seq_data_list]
        all_seq_bytes = b''.join(seq_data_list)
        all_seq_array = np.frombuffer(all_seq_bytes, dtype=np.uint8)
        
        # Apply base mapping to ENTIRE batch at once
        if keep_bases:
            mapped_batch = all_seq_array
        elif binary_map is not None:
            mapped_batch = binary_map[all_seq_array]
        else:
            mapped_batch = base_map[all_seq_array]
        
        cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
        
        # Create views (instead of copies for major speedup)
        seq_arrays = [mapped_batch[cumsum_lengths[i]:cumsum_lengths[i+1]] for i in range(batch_size)]
        
        # Process qualities (if needed)
        qual_arrays = []
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
            
            # Fast path: skip pair detection if not needed
            pair_number = 0
            if paired_end:
                # Inline detection (avoids function call)
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