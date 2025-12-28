import numpy as np
from typing import List, Optional
from data_structures import FASTQRecord
from quality_processing import apply_quality_to_bases


# Reserved sequence start marker as val of 255
SEQUENCE_START_MARKER = np.uint8(255)

def process_and_write_records(records: List[FASTQRecord], outfile, base_map: np.ndarray,
                               quality_scaling: str, custom_formula: Optional[str],
                               phred_alphabet_max: int, min_quality: int, keep_bases: bool,
                               binary: bool, keep_quality: bool,
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
        
        if not keep_bases:
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

        if not keep_bases:
            outfile.write(bytes([SEQUENCE_START_MARKER])) # Reserved 255 val

        seq_data = processed_seqs[idx]
        
        if not keep_bases and np.any(seq_data == 255): # Don't want to write out 255
            # This should never happen with proper base mapping, but safety check yk
            print(f"WARNING: Sequence {idx} contains value 255, clamping to 254")
            seq_data = np.clip(seq_data, 0, 254)

        # Write sequence data in appropriate format
        if binary:
            outfile.write(seq_data.tobytes())
        else:
            if keep_bases:
                outfile.write(seq_data.tobytes())  # Write ASCII directly
            else:
                outfile.write(b''.join(BYTE_LOOKUP[seq_data]))  # Convert numbers to text
        
        outfile.write(b'\n')
        
        # Write quality scores if keep_quality is enabled
        if keep_quality and len(record.quality) > 0:
            outfile.write(b'+\n')
            outfile.write(record.quality_string)  # ASCII quality string
            outfile.write(b'\n')