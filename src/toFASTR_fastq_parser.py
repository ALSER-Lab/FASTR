from typing import Dict, List, Optional, Tuple

import numpy as np

from data_structures import FASTQRecord
from toFASTR_header_compression import (adaptive_compress_header,
                                        analyze_headers_for_pattern,
                                        compress_header)


def parse_fastq_records_from_buffer(
    buffer: bytes,
    start_index: int,
    base_map: np.ndarray,
    phred_map: Optional[np.ndarray],
    compress_headers: bool,
    sequencer_type: str,
    keep_bases: bool,
    keep_quality: bool,
    adaptive_structure: Optional[str] = None,
    adaptive_delimiter: Optional[str] = None,
    adaptive_sample_size: int = 100,
) -> Tuple[List[FASTQRecord], bytes, Optional[Dict], Optional[str], Optional[str], int]:
    """
    Parse FASTQ records from a buffer using vectorized operations.
    Returns (records, leftover_buffer, current_flowcell_metadata, structure_template, delimiter, num_records)
    """

    records = []
    current_flowcell_metadata = None
    structure_template = adaptive_structure
    delimiter = adaptive_delimiter
    lines = buffer.split(b"\n")

    if lines and lines[-1] == b"":
        lines = lines[:-1]

    record_boundaries = []
    i = 0
    while i < len(lines):
        if lines[i].startswith(b"@"):
            header_idx = i
            i += 1
            seq_start_idx = i

            while i < len(lines) and not lines[i].startswith(b"+"):
                i += 1

            if i >= len(lines):
                break  # Incomplete record

            plus_idx = i
            i += 1
            qual_start_idx = i
            seq_length = sum(len(lines[j]) for j in range(seq_start_idx, plus_idx))
            qual_length = 0
            while i < len(lines) and qual_length < seq_length:
                qual_length += len(lines[i])
                i += 1

            if qual_length >= seq_length:
                record_boundaries.append(
                    (header_idx, seq_start_idx, plus_idx, qual_start_idx, i)
                )
            else:
                break  # Incomplete record
        else:
            i += 1

    num_complete_records = len(record_boundaries)

    if num_complete_records == 0:
        return (
            records,
            buffer,
            current_flowcell_metadata,
            structure_template,
            delimiter,
            0,
        )

    leftover_start = record_boundaries[-1][-1] if record_boundaries else 0
    leftover = b"\n".join(lines[leftover_start:])

    # If using adaptive mode and we don't have a structure yet, analyze headers first
    if compress_headers and sequencer_type == "adaptive" and structure_template is None:
        # Extract sample of headers for analysis
        sample_headers = []
        max_samples = min(
            num_complete_records * 4, adaptive_sample_size * 4
        )  # User-inputted sampling size
        for i in range(
            0, max_samples, 4
        ):  # Divide by 4 bc each "block" of fastq has 4 lines
            header_str = lines[i].decode("utf-8", errors="ignore")
            if header_str.startswith("@"):
                header_str = header_str[1:]  # Remove @ prefix
            sample_headers.append(header_str)

        structure_template, delimiter, current_flowcell_metadata = (
            analyze_headers_for_pattern(
                sample_headers, sample_size=adaptive_sample_size
            )
        )

    # Pre-create binary map ONCE if needed

    # Process in batches to avoid numpy split overhead
    BATCH_SIZE = 10000  # Process 10k records at a time

    for batch_start in range(0, num_complete_records, BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, num_complete_records)
        batch_size = batch_end - batch_start

        header_strs = []
        seq_data_list = []
        qual_data_list = []

        for rec_idx in range(batch_start, batch_end):
            header_idx, seq_start_idx, plus_idx, qual_start_idx, qual_end_idx = (
                record_boundaries[rec_idx]
            )

            header_strs.append(lines[header_idx].decode("utf-8", errors="ignore"))

            seq_lines = [
                lines[j].replace(b"\r", b"") for j in range(seq_start_idx, plus_idx)
            ]  # Join sequences that may be multi line
            seq_data_list.append(b"".join(seq_lines))

            qual_lines = [
                lines[j].replace(b"\r", b"")
                for j in range(qual_start_idx, qual_end_idx)
            ]  # Join quality that may be multi line
            qual_data_list.append(b"".join(qual_lines))

        seq_lengths = [len(s) for s in seq_data_list]
        all_seq_bytes = b"".join(seq_data_list)
        all_seq_array = np.frombuffer(all_seq_bytes, dtype=np.uint8)

        # Apply base mapping to ENTIRE batch at once
        if keep_bases:
            mapped_batch = all_seq_array
        else:
            mapped_batch = base_map[all_seq_array]

        cumsum_lengths = np.concatenate([[0], np.cumsum(seq_lengths)])

        # Create views (instead of copies for major speedup)
        seq_arrays = [
            mapped_batch[cumsum_lengths[i] : cumsum_lengths[i + 1]]
            for i in range(batch_size)
        ]

        # Process qualities (if needed)
        qual_arrays = []
        if phred_map is not None:
            all_qual_bytes = b"".join(qual_data_list)
            all_qual_array = np.frombuffer(all_qual_bytes, dtype=np.uint8)
            qual_mapped = phred_map[all_qual_array]
            # Use views for quality too
            qual_arrays = [
                qual_mapped[cumsum_lengths[i] : cumsum_lengths[i + 1]]
                for i in range(batch_size)
            ]

        # Build records
        for idx in range(batch_size):
            record_index = start_index + batch_start + idx
            header_str = header_strs[idx]

            header_idx_in_rec = record_boundaries[batch_start + idx][0]
            original_header = lines[header_idx_in_rec] + b"\n"

            # Header compression
            if compress_headers and sequencer_type != "none":
                if sequencer_type == "adaptive":
                    # Use adaptive compression with detected structure
                    common, unique_id = adaptive_compress_header(
                        header_str,
                        structure_template,
                        delimiter,
                        current_flowcell_metadata,
                    )
                else:
                    # Use standard sequencer based compression
                    common, unique_id, structure = compress_header(
                        header_str, sequencer_type
                    )
                    if common:
                        current_flowcell_metadata = common

                    if structure and structure_template is None:
                        structure_template = structure

                header = f"@{unique_id}\n".encode("utf-8")
            else:
                header = original_header

            # Get quality
            qual_vals = (
                qual_arrays[idx] if qual_arrays else np.array([], dtype=np.uint8)
            )

            record = FASTQRecord(
                index=record_index,
                header=header,
                sequence=seq_arrays[idx],
                quality=qual_vals,
                quality_string=qual_data_list[idx] if keep_quality else b"",
                original_header=original_header,
            )
            records.append(record)

    return (
        records,
        leftover,
        current_flowcell_metadata,
        structure_template,
        delimiter,
        num_complete_records,
    )

