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


def process_and_write_records(
    records: List[FASTQRecord],
    outfile,
    base_map: np.ndarray,
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
    quality_lookup_table=None,
):
    """
    Apply quality transformations to records and write to FASTR.
    """
    if not records:
        return

    seq_lengths = [len(r.sequence) for r in records]
    flat_mapped = np.concatenate([r.sequence for r in records]).astype(np.uint8)
    all_quals = [r.quality for r in records if len(r.quality) > 0]

    if len(all_quals) > 0 and not keep_bases:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        flat_mapped = apply_quality_to_bases(
            flat_mapped,
            flat_quals,
            base_map,
            quality_lookup_table=quality_lookup_table,
        )

    if min_quality > 0 and len(all_quals) > 0:
        flat_quals = np.concatenate(all_quals).astype(np.uint8)
        low_quality_mask = flat_quals < min_quality
        if keep_bases:
            flat_mapped[low_quality_mask] = ord("N")
        else:
            flat_mapped[low_quality_mask] = base_map[ord("N")]

    cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])
    processed_seqs = [
        flat_mapped[cumsum_lengths[i] : cumsum_lengths[i + 1]]
        for i in range(len(records))
    ]

    for idx, record in enumerate(records):
        seq_data = processed_seqs[idx]
        seq_data[seq_data == 10] = 255

        if headers_buffer is not None:
            headers_buffer.write(record.original_header)

        if mode in [0, 1, 2]:
            if not remove_repeating_header:
                header_bytes = record.header
                if header_bytes.endswith(b"\n"):
                    header_bytes = header_bytes[:-1]
                outfile.write(header_bytes)
                outfile.write(b"\n")

        if binary:
            outfile.write(seq_data.tobytes())
            outfile.write(b"\n")
        else:
            if keep_bases:
                outfile.write(seq_data.tobytes())
            else:
                outfile.write(b"".join(BYTE_LOOKUP[seq_data]))
            outfile.write(b"\n")

        if keep_quality and len(record.quality) > 0:
            outfile.write(b"+\n")
            outfile.write(record.quality_string)
            outfile.write(b"\n")
