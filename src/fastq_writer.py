import numpy as np
from typing import List, Optional
from data_structures import FASTQRecord
from quality_processing import apply_quality_to_bases


def process_and_write_records(records: List[FASTQRecord], outfile, base_map: np.ndarray,
                               quality_scaling: str, custom_formula: Optional[str],
                               phred_alphabet_max: int, min_quality: int, keep_bases: bool,
                               binary_bases: bool, binary: bool, keep_quality: bool,
                               remove_repeating_header: bool, compress_headers: bool,
                               BYTE_LOOKUP: np.ndarray):
    """
    Process quality scaling on records and write to file.
    Applies quality-based transformations and outputs in specified format.
    """
    if not records:
        return
    
    # Store sequence lengths before concatenation
    seq_lengths = [len(r.sequence) for r in records]
    flat_mapped = np.concatenate([r.sequence for r in records]).astype(np.uint8)
    
    # Collect all quality scores
    all_quals = [r.quality for r in records if len(r.quality) > 0]
    
    # Apply quality scaling if needed
    if len(all_quals) > 0 and quality_scaling != 'none':
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        
        if not keep_bases and not binary_bases:
            # Apply scaling on flattened arrays
            flat_mapped = apply_quality_to_bases(
                flat_mapped, flat_quals, base_map, quality_scaling,
                custom_formula, phred_alphabet_max
            )
    
    # Apply minimum quality filtering
    if min_quality > 0 and len(all_quals) > 0:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        low_quality_mask = flat_quals < min_quality
        
        if keep_bases:
            flat_mapped[low_quality_mask] = ord('N')  # Keep as ASCII
        elif binary_bases:
            flat_mapped[low_quality_mask] = 4  # N = 4 in binary encoding
        else:
            flat_mapped[low_quality_mask] = base_map[ord('N')]  # Convert to numeric
    
    # Split back into individual sequences
    cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
    processed_seqs = [
        flat_mapped[cumsum_lengths[i]:cumsum_lengths[i+1]]
        for i in range(len(records))
    ]
    
    # Write all sequences
    for idx, record in enumerate(records):
        # Only write headers if remove_repeating_header is disabled
        if not remove_repeating_header:
            outfile.write(record.header)
        
        # Add start marker for sequence
        outfile.write(b'<')
        
        seq_data = processed_seqs[idx]
        
        # Write sequence data in appropriate format
        if binary:
            outfile.write(seq_data.tobytes())
        else:
            if keep_bases:
                outfile.write(seq_data.tobytes())  # Write ASCII directly
            elif binary_bases:
                outfile.write(b''.join(BYTE_LOOKUP[seq_data]))  # Write binary bases as text
            else:
                outfile.write(b''.join(BYTE_LOOKUP[seq_data]))  # Convert numbers to text
        
        outfile.write(b'\n')
        
        # Write quality scores if keep_quality is enabled
        if keep_quality and len(record.quality) > 0:
            outfile.write(b'+\n')
            outfile.write(record.quality)  # ASCII quality string
            outfile.write(b'\n')