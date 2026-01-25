import logging
from typing import List, Optional

import numpy as np

from data_structures import FASTQRecord
from toFASTR_quality_processing import apply_quality_to_bases

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Reserved sequence start marker as val of 255
SEQUENCE_START_MARKER = np.uint8(255)


def process_and_write_records(
    records: List[FASTQRecord],
    outfile,
    base_map: np.ndarray,
    quality_scaling: str,
    custom_formula: Optional[str],
    phred_alphabet_max: int,
    min_quality: int,
    keep_bases: bool,
    binary: bool,
    keep_quality: bool,
    remove_repeating_header: bool,
    compress_headers: bool,
    BYTE_LOOKUP: np.ndarray,
    headers_buffer=None,
    mode: int = None,
    safe_mode: bool = False,
):
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
    if len(all_quals) > 0 and quality_scaling != "none":
        flat_quals = np.concatenate(all_quals).astype(np.uint8)

        if not keep_bases:
            # Apply scaling on flattened arrays
            flat_mapped = apply_quality_to_bases(
                flat_mapped,
                flat_quals,
                base_map,
                quality_scaling,
                custom_formula,
                phred_alphabet_max,
            )

    # Apply minimum quality filtering
    if min_quality > 0 and len(all_quals) > 0:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        low_quality_mask = flat_quals < min_quality

        if keep_bases:
            flat_mapped[low_quality_mask] = ord("N")  # Keep as ASCII
        else:
            flat_mapped[low_quality_mask] = base_map[ord("N")]  # Convert to numeric

    # Split back into individual sequences
    cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
    processed_seqs = [
        flat_mapped[cumsum_lengths[i] : cumsum_lengths[i + 1]]
        for i in range(len(records))
    ]

    # Write all sequences
    for idx, record in enumerate(records):
        seq_data = processed_seqs[idx]

        if not keep_bases and np.any(seq_data >= 255):
            logger.warning(
                f"Sequence {idx} contains reserved value 255 (likely from custom base ranges). "
                f"Clamping to 254. Check your --gray_* arguments to avoid data loss."
            )
            seq_data = np.clip(seq_data, 0, 254)

        # Write headers to separate file if headers_buffer provided (Mode 3)
        if headers_buffer is not None:
            headers_buffer.write(record.original_header)
            # For mode 3 we write 255 before bases
            if binary:
                outfile.write(bytes([SEQUENCE_START_MARKER]))

        if mode == 0:
            if not remove_repeating_header:
                header_bytes = record.header
                if header_bytes.endswith(b"\n"):
                    header_bytes = header_bytes[:-1]
                outfile.write(header_bytes)
                outfile.write(b"\n")

        # Mode 2 and mode 1 we write 255 after header
        elif mode in [1, 2]:
            if not remove_repeating_header:
                header_bytes = record.header
                if header_bytes.endswith(b"\n"):
                    header_bytes = header_bytes[:-1]
                outfile.write(header_bytes)
                if safe_mode:  # 255 then \n
                    outfile.write(bytes([SEQUENCE_START_MARKER]))
                    outfile.write(b"\n")
                else:  # Just \n
                    outfile.write(b"\n")

        # Write sequence data
        if binary:
            outfile.write(seq_data.tobytes())
            # Mode 3 doesn't have newlines between sequences, other modes do
            outfile.write(b"\n")
        else:
            if keep_bases:
                outfile.write(seq_data.tobytes())  # Write ASCII directly
            else:
                outfile.write(
                    b"".join(BYTE_LOOKUP[seq_data])
                )  # Convert numbers to text
            outfile.write(b"\n")

        # Write quality scores if keep_quality is enabled
        if keep_quality and len(record.quality) > 0:
            outfile.write(b"+\n")
            outfile.write(record.quality_string)
            outfile.write(b"\n")

