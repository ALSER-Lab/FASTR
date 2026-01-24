import io
import logging
import os
import traceback
from typing import List, Optional, Tuple

import numpy as np

from base_mapping import create_base_map, reverse_base_map
from data_structures import MetadataBlock
from header_reconstruction import reconstruct_header_from_structure
from metadata_parser import find_structure_prefix
from quality_reconstruction import quality_to_ascii, reverse_scaling_to_quality


def process_chunk_worker_reconstruction(
    chunk_data,
    mmap_path,
    reverse_map,
    subtract_table,
    metadata_blocks,
    inverse_tables,
    phred_alphabet_max,
    phred_offset,
    sra_acc,
    mode,
    length_flag=False,
    headers_file_path=None,
    headers_mmap_info=None,
    second_head_flag=False,
    safe_mode_flag=False,
    data_start_byte=None,
):
    """
    Worker function for parallel processing of FASTR chunks.
    Handles modes 0-3: (0) headers only, (1) bases only, (2) full reconstruction, (3) no repeating headers.
    Returns: (chunk_id, output_bytes, sequence_count)

    """
    try:
        chunk_id, abs_start, abs_end, start_seq_idx = chunk_data
        logger = logging.getLogger(__name__)
        if mode == 3 and headers_mmap_info:
            if not hasattr(process_chunk_worker_reconstruction, "_header_index"):
                mmap_path_hdr, shape, dtype = headers_mmap_info
                process_chunk_worker_reconstruction._header_index = np.memmap(
                    mmap_path_hdr, dtype=dtype, mode="r", shape=shape
                )
                logger.debug(f"Worker loaded header index: {shape[0]:,} offsets")
            header_index = process_chunk_worker_reconstruction._header_index
        else:
            header_index = None

        MAX_CHUNK_EXTENSION = 10 * 1024 * 1024

        with open(mmap_path, "rb") as f:
            f.seek(abs_start)
            chunk_read_size = (abs_end - abs_start) + MAX_CHUNK_EXTENSION
            data = f.read(chunk_read_size)

        chunk_size = abs_end - abs_start

        if mode == 0:
            output_bytes, sequences_processed = _process_mode_0(
                data=data,
                chunk_size=chunk_size,
                start_seq_idx=start_seq_idx,
                reverse_map=reverse_map,
                subtract_table=subtract_table,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                length_flag=length_flag,
                second_head_flag=second_head_flag,
            )
        elif mode == 1:
            output_bytes, sequences_processed = _process_mode_1(
                data=data,
                chunk_size=chunk_size,
                start_seq_idx=start_seq_idx,
                reverse_map=reverse_map,
                subtract_table=subtract_table,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                second_head_flag=second_head_flag,
                safe_mode_flag=safe_mode_flag,
            )
        elif mode == 2:
            output_bytes, sequences_processed = _process_mode_2(
                data=data,
                chunk_size=chunk_size,
                start_seq_idx=start_seq_idx,
                reverse_map=reverse_map,
                subtract_table=subtract_table,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                sra_accession=sra_acc,
                length_flag=length_flag,
                second_head_flag=second_head_flag,
                safe_mode_flag=safe_mode_flag,
            )
        elif mode == 3:
            output_bytes, sequences_processed = _process_mode_3(
                data=data,
                chunk_size=chunk_size,
                start_seq_idx=start_seq_idx,
                reverse_map=reverse_map,
                subtract_table=subtract_table,
                metadata_blocks=metadata_blocks,
                inverse_tables=inverse_tables,
                phred_alphabet_max=phred_alphabet_max,
                phred_offset=phred_offset,
                sra_accession=sra_acc,
                headers_file_path=headers_file_path,
                header_index=header_index,
                second_head_flag=second_head_flag,
            )
        else:
            raise ValueError(f"Invalid mode: {mode}")

        return (chunk_id, output_bytes, sequences_processed)

    except Exception as e:
        logger.warning(f"ERROR in worker processing chunk {chunk_id}: {e}")
        traceback.print_exc()
        raise


def _process_mode_0(
    data: bytes,
    chunk_size: int,
    start_seq_idx: int,
    reverse_map: np.ndarray,
    subtract_table: np.ndarray,
    metadata_blocks: List[MetadataBlock],
    inverse_tables: List[np.ndarray],
    phred_alphabet_max: int,
    phred_offset: int,
    length_flag: bool,
    second_head_flag: bool,
) -> Tuple[bytes, int]:
    """
    Process Mode 0: Headers compressed, bases and quality as plain text.
    Reconstructs full headers from compressed unique IDs using structure templates.

    Returns: (output_bytes, sequences_processed)
    """
    output_buffer = io.BytesIO()
    sequence_count = start_seq_idx
    current_metadata_idx = 0
    current_metadata = metadata_blocks[0] if metadata_blocks else None
    MAX_CHUNK_EXTENSION = 10 * 1024 * 1024
    sequences_in_chunk = 0
    cursor = 0

    while cursor < chunk_size:
        if data[cursor : cursor + 1] != b"@":
            # Try to find next @
            next_at = data.find(
                b"@", cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION)
            )
            if next_at == -1:
                break
            cursor = next_at

        line1_end = data.find(b"\n", cursor)  # Header
        if line1_end == -1 or line1_end >= len(data):
            break

        compressed_header = (
            data[cursor + 1 : line1_end].decode("utf-8", errors="replace").strip()
        )
        pair_number = 0

        line2_start = line1_end + 1
        seq_cursor = line2_start
        bases_len = 0
        line3_start = -1
        bases_parts = []

        while seq_cursor < len(data):
            line_end = data.find(b"\n", seq_cursor)
            if line_end == -1:
                break

            next_line_start = line_end + 1
            if (
                next_line_start < len(data)
                and data[next_line_start : next_line_start + 1] == b"+"
            ):
                bases_parts.append(data[seq_cursor:line_end])
                bases_len = sum(len(p) for p in bases_parts)
                line3_start = next_line_start
                break

            bases_parts.append(data[seq_cursor:line_end])
            seq_cursor = line_end + 1

        if line3_start == -1:
            break

        bases_line = b"".join(bases_parts).decode("utf-8", errors="replace").strip()
        line3_end = data.find(b"\n", line3_start)  # Line 3
        if line3_end == -1:
            break

        plus_line = (
            data[line3_start:line3_end].decode("utf-8", errors="replace").strip()
        )

        if not plus_line.startswith("+"):
            cursor = line3_end + 1
            continue

        line4_start = line3_end + 1
        line4_end = line4_start + bases_len

        if line4_end > len(data):
            break

        quality_line = data[line4_start:line4_end].decode("utf-8", errors="replace")

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
                pair_number,
            )
        else:
            full_header = f"@{unique_id}"
            if pair_number > 0:
                full_header = f"{full_header}/{pair_number}"

        # Add length flag if needed
        if length_flag:
            full_header = f"{full_header} length={len(bases_line)}"

        # Write reconstructed record
        output_buffer.write(full_header.encode("utf-8"))
        output_buffer.write(b"\n")
        output_buffer.write(bases_line.encode("utf-8"))
        output_buffer.write(b"\n")

        if second_head_flag:
            output_buffer.write(b"+")
            output_buffer.write(full_header[1:].encode("utf-8"))  # Remove '@'
            output_buffer.write(b"\n")
        else:
            output_buffer.write(b"+\n")

        output_buffer.write(quality_line.encode("utf-8"))
        output_buffer.write(b"\n")

        sequence_count += 1
        sequences_in_chunk += 1

        # Move cursor past this record (skip trailing newline if present)
        cursor = line4_end
        if cursor < len(data) and data[cursor : cursor + 1] == b"\n":
            cursor += 1

    return (output_buffer.getvalue(), sequences_in_chunk)


def _process_mode_1(
    data: bytes,
    chunk_size: int,
    start_seq_idx: int,
    reverse_map: np.ndarray,
    subtract_table: np.ndarray,
    metadata_blocks: List[MetadataBlock],
    inverse_tables: List[np.ndarray],
    phred_alphabet_max: int,
    phred_offset: int,
    second_head_flag: bool,
    safe_mode_flag: bool,
) -> Tuple[bytes, int]:
    """
       Process Mode 1: Only bases compressed (original headers intact).
       Decodes compressed bases back to ACGT and reconstructs quality scores.
       Handles both safe mode (with \\xff validation) and unsafe mode parsing.

    Returns: (output_bytes, sequences_processed)
    """
    output_buffer = io.BytesIO()
    sequence_count = start_seq_idx
    current_metadata_idx = 0
    current_metadata = metadata_blocks[0] if metadata_blocks else None
    current_inverse_table = inverse_tables[0] if inverse_tables else None
    MAX_CHUNK_EXTENSION = 10 * 1024 * 1024
    sequences_in_chunk = 0
    cursor = 0

    while cursor < chunk_size:
        if safe_mode_flag:
            # Safe mode, format is: \n@header(255)\nbases\n
            # Find next @ with validation
            if cursor == 0 and data[0:1] == b"@":
                # First sequence in chunk starts at position 0
                header_start_rel = 0
            else:
                # Always look for \n@ for subsequent sequences
                header_start_rel = data.find(b"\n@", cursor)
                if header_start_rel != -1:
                    header_start_rel += 1

            if header_start_rel < 0 or header_start_rel >= chunk_size:
                break

            candidate = header_start_rel  # Then validate by finding \xff after this

            while candidate >= 0:
                next_xff_rel = data.find(
                    b"\xff",
                    candidate,
                    min(len(data), chunk_size + MAX_CHUNK_EXTENSION),
                )
                next_at_rel = data.find(
                    b"\n@",
                    candidate + 1,
                    min(len(data), chunk_size + MAX_CHUNK_EXTENSION),
                )
                # If another \n@ comes before \xff, current candidate is invalid
                if next_at_rel != -1 and (
                    next_xff_rel == -1 or next_at_rel < next_xff_rel
                ):
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

            header = (
                data[header_start:xff_pos].decode("utf-8", errors="replace").strip()
            )

            if xff_pos + 1 < len(data) and data[xff_pos + 1 : xff_pos + 2] == b"\n":
                seq_start_rel = xff_pos + 2
            else:
                seq_start_rel = xff_pos + 1

            candidate_end_rel = seq_start_rel
            seq_end_rel = -1

            while True:
                next_at_rel = data.find(
                    b"\n@",
                    candidate_end_rel,
                    min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION),
                )

                if next_at_rel != -1:
                    next_xff_rel = data.find(
                        b"\xff",
                        next_at_rel + 1,
                        min(len(data), next_at_rel + 1000),
                    )
                    following_at_rel = data.find(
                        b"\n@",
                        next_at_rel + 2,
                        min(len(data), next_at_rel + 1000),
                    )

                    if following_at_rel != -1 and (
                        next_xff_rel == -1 or following_at_rel < next_xff_rel
                    ):
                        candidate_end_rel = following_at_rel + 1
                        if candidate_end_rel > seq_start_rel + MAX_CHUNK_EXTENSION:
                            break
                        continue

                    seq_end_rel = next_at_rel
                    break
                else:
                    seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                    while seq_end_rel > seq_start_rel and data[
                        seq_end_rel - 1 : seq_end_rel
                    ] in (b"\n", b"\r", b" "):
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
            header_start_rel = data.find(
                b"@", cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION)
            )
            if header_start_rel == -1:
                break

            header_end_rel = data.find(
                b"\n",
                header_start_rel,
                min(len(data), header_start_rel + 10000),
            )
            if header_end_rel == -1:
                break

            header = data[header_start_rel:header_end_rel].decode(
                "utf-8", errors="ignore"
            )

            seq_start_rel = header_end_rel + 1

            next_at_rel = data.find(
                b"@",
                seq_start_rel + 1,
                min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION),
            )

            if next_at_rel != -1 and next_at_rel > seq_start_rel:
                if data[next_at_rel - 1 : next_at_rel] == b"\n":
                    seq_end_rel = next_at_rel - 1
                else:
                    seq_end_rel = next_at_rel
            else:
                seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                if (
                    seq_end_rel > seq_start_rel
                    and data[seq_start_rel:seq_end_rel][-1:] == b"\n"
                ):
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
        bases = bases_array.tobytes().decode("ascii")

        # Reconstruct quality
        if current_inverse_table is not None:
            quality_scores = reverse_scaling_to_quality(
                seq_array,
                subtract_table,
                current_inverse_table,
                max_phred=phred_alphabet_max,
            )
            quality_string = quality_to_ascii(quality_scores, phred_offset).decode(
                "ascii"
            )

        # Write output
        output_header = (
            f"{header}\n{bases}\n+\n"
            if not second_head_flag
            else f"{header}\n{bases}\n+{header[1:]}\n"
        )
        output_buffer.write(output_header.encode("ascii"))
        output_buffer.write(quality_string.encode("ascii"))
        output_buffer.write(b"\n")
        sequence_count += 1
        sequences_in_chunk += 1

    return (output_buffer.getvalue(), sequences_in_chunk)


def _process_mode_2(
    data: bytes,
    chunk_size: int,
    start_seq_idx: int,
    reverse_map: np.ndarray,
    subtract_table: np.ndarray,
    metadata_blocks: List[MetadataBlock],
    inverse_tables: List[np.ndarray],
    phred_alphabet_max: int,
    phred_offset: int,
    sra_accession: Optional[str],
    length_flag: bool,
    second_head_flag: bool,
    safe_mode_flag: bool,
) -> Tuple[bytes, int]:
    """
    Process Mode 2: Full reconstruction (headers + bases compressed).
    Parses compressed unique_id from headers, extracts pair numbers,
    reconstructs full headers from templates, and decodes bases/quality.

    Returns: (output_bytes, sequences_processed)
    """
    output_buffer = io.BytesIO()
    sequence_count = start_seq_idx
    current_metadata_idx = 0
    current_metadata = metadata_blocks[0] if metadata_blocks else None
    current_inverse_table = inverse_tables[0] if inverse_tables else None
    MAX_CHUNK_EXTENSION = 10 * 1024 * 1024
    # Safe mode, format is: \n@header(255)\nbases\n
    sequences_in_chunk = 0
    cursor = 0

    while cursor < chunk_size:
        if safe_mode_flag:
            # Find next @ with validation
            if cursor == 0 and data[0:1] == b"@":
                # First sequence in chunk starts at position 0
                header_start_rel = 0
            else:
                # Always look for \n@ for subsequent sequences
                header_start_rel = data.find(b"\n@", cursor)
                if header_start_rel != -1:
                    header_start_rel += 1  # Skip the \n, point to the @

            if header_start_rel < 0 or header_start_rel >= chunk_size:
                break

            candidate = header_start_rel  # Then validate by finding \xff after this

            while candidate >= 0:
                next_xff_rel = data.find(
                    b"\xff",
                    candidate,
                    min(len(data), chunk_size + MAX_CHUNK_EXTENSION),
                )
                next_at_rel = data.find(
                    b"\n@",
                    candidate + 1,
                    min(len(data), chunk_size + MAX_CHUNK_EXTENSION),
                )

                # If another \n@ comes before \xff, current candidate is invalid
                if next_at_rel != -1 and (
                    next_xff_rel == -1 or next_at_rel < next_xff_rel
                ):
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
            header_content = (
                data[header_start + 1 : xff_pos]
                .decode("utf-8", errors="replace")
                .strip()
            )

            # Bases start after \xff\n or \xff
            if xff_pos + 1 < len(data) and data[xff_pos + 1 : xff_pos + 2] == b"\n":
                seq_start_rel = xff_pos + 2
            else:
                seq_start_rel = xff_pos + 1

            # Find sequence end (next valid @ header)
            candidate_end_rel = seq_start_rel
            seq_end_rel = -1

            while True:
                next_at_rel = data.find(
                    b"\n@",
                    candidate_end_rel,
                    min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION),
                )

                if next_at_rel != -1:
                    # Validate our \n@
                    next_xff_rel = data.find(
                        b"\xff",
                        next_at_rel + 1,
                        min(len(data), next_at_rel + 1000),
                    )
                    following_at_rel = data.find(
                        b"\n@",
                        next_at_rel + 2,
                        min(len(data), next_at_rel + 1000),
                    )

                    if following_at_rel != -1 and (
                        next_xff_rel == -1 or following_at_rel < next_xff_rel
                    ):
                        candidate_end_rel = following_at_rel + 1
                        if candidate_end_rel > seq_start_rel + MAX_CHUNK_EXTENSION:
                            break
                        continue

                    # Valid header if sequence ends at the \n
                    seq_end_rel = next_at_rel
                    break
                else:
                    seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                    while seq_end_rel > seq_start_rel and data[
                        seq_end_rel - 1 : seq_end_rel
                    ] in (b"\n", b"\r", b" "):
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
            start_idx_rel = data.find(
                b"@", cursor, min(len(data), chunk_size + MAX_CHUNK_EXTENSION)
            )
            if start_idx_rel == -1:
                break

            header_end_rel = data.find(
                b"\n", start_idx_rel, min(len(data), start_idx_rel + 10000)
            )
            if header_end_rel == -1:
                break

            header_content = (
                data[start_idx_rel + 1 : header_end_rel]
                .decode("utf-8", errors="ignore")
                .strip()
            )

            seq_start_rel = header_end_rel + 1

            next_header_rel = data.find(
                b"\n@",
                seq_start_rel,
                min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION),
            )

            if next_header_rel != -1:
                seq_end_rel = next_header_rel
            else:
                seq_end_rel = min(len(data), seq_start_rel + MAX_CHUNK_EXTENSION)
                while seq_end_rel > seq_start_rel and data[
                    seq_end_rel - 1 : seq_end_rel
                ] in (b"\n", b"\r", b" ", b"\t"):
                    seq_end_rel -= 1

            seq_data = data[seq_start_rel:seq_end_rel]
            cursor = seq_end_rel + 1

        unique_id = None
        pair_number = 0

        if "/" in header_content:
            parts = header_content.rsplit("/", 1)
            if (
                len(parts) == 2
                and parts[1].isdigit()
                and len(parts[1]) == 1
                and int(parts[1]) in [1, 2]
            ):
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
                pair_number,
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
        bases = bases_array.tobytes().decode("ascii")

        if current_inverse_table is not None:
            quality_scores = reverse_scaling_to_quality(
                seq_array,
                subtract_table,
                current_inverse_table,
                max_phred=phred_alphabet_max,
            )
            quality_string = quality_to_ascii(quality_scores, phred_offset).decode(
                "ascii"
            )

        if length_flag:
            header = f"{header} length={len(bases)}"

        output_buffer.write(header.encode("utf-8"))
        output_buffer.write(b"\n")
        output_buffer.write(bases.encode("ascii"))
        output_buffer.write(b"\n+")
        if second_head_flag:
            output_buffer.write(header[1:].encode("utf-8"))
        output_buffer.write(b"\n")
        output_buffer.write(quality_string.encode("ascii"))
        output_buffer.write(b"\n")

        sequence_count += 1
        sequences_in_chunk += 1

    return (output_buffer.getvalue(), sequences_in_chunk)


def _process_mode_3(
    data: bytes,
    chunk_size: int,
    start_seq_idx: int,
    reverse_map: np.ndarray,
    subtract_table: np.ndarray,
    metadata_blocks: List[MetadataBlock],
    inverse_tables: List[np.ndarray],
    phred_alphabet_max: int,
    phred_offset: int,
    sra_accession: Optional[str],
    headers_file_path: Optional[str],
    header_index: Optional[np.ndarray],
    second_head_flag: bool,
) -> Tuple[bytes, int]:
    """
    Process Mode 3: No repeating headers (requires external headers file).
    Uses \\xff markers as sequence delimiters, loads headers from external file
    with batch caching for performance. Falls back to generated headers on cache miss.

    Returns: (output_bytes, sequences_processed)
    """
    output_buffer = io.BytesIO()
    sequence_count = start_seq_idx
    current_metadata_idx = 0
    current_metadata = metadata_blocks[0] if metadata_blocks else None
    current_inverse_table = inverse_tables[0] if inverse_tables else None
    logger = logging.getLogger(__name__)
    sequences_in_chunk = 0
    MAX_CHUNK_EXTENSION = 10 * 1024 * 1024

    headers_cache = {}
    if header_index is not None and headers_file_path:
        search_region_size = min(chunk_size + 1000000, len(data))
        actual_marker_count = data[:search_region_size].count(
            b"\xff"
        )  # Count how many sequences in our chunk based on delimiter

        batch_start_idx = start_seq_idx
        batch_end_idx = min(start_seq_idx + actual_marker_count, len(header_index))

        if batch_start_idx < len(header_index):
            start_offset = int(header_index[batch_start_idx])

            if batch_end_idx < len(header_index):
                end_offset = int(header_index[batch_end_idx])
            else:
                end_offset = os.path.getsize(headers_file_path)

            try:
                with open(headers_file_path, "rb") as hf:
                    hf.seek(start_offset)
                    batch_data = hf.read(end_offset - start_offset)

                header_lines = (
                    batch_data.decode("utf-8", errors="ignore").strip().split("\n")
                )
                for i, header_line in enumerate(header_lines):
                    cache_idx = batch_start_idx + i
                    if cache_idx < len(header_index):
                        headers_cache[cache_idx] = header_line.strip()

            except Exception as e:
                logger.debug(f"Header load failed: {e}")

    search_limit = min(len(data), chunk_size + MAX_CHUNK_EXTENSION)
    data_array = np.frombuffer(data[:search_limit], dtype=np.uint8)
    all_marker_positions = np.where(data_array == 255)[0]

    marker_positions = all_marker_positions[all_marker_positions < chunk_size]

    if len(marker_positions) == 0:
        return (output_buffer.getvalue(), 0)
    for i in range(len(marker_positions)):
        marker_pos = int(marker_positions[i])
        seq_start_rel = marker_pos + 1

        next_marker_idx = i + 1
        if next_marker_idx < len(all_marker_positions):
            seq_end_rel = int(all_marker_positions[next_marker_idx])
        else:
            seq_end_rel = len(data)

        if (
            seq_end_rel > seq_start_rel
            and data[seq_end_rel - 1 : seq_end_rel]
            == b"\n"  # Exclude the last character (always a \n from our conversions)
        ):
            seq_end_rel -= 1

        if seq_end_rel <= seq_start_rel:
            continue

        seq_data = data[seq_start_rel:seq_end_rel]

        if len(seq_data) == 0:
            continue

        # Check for metadata block changes
        if len(metadata_blocks) > 1 and current_metadata_idx < len(metadata_blocks) - 1:
            next_metadata = metadata_blocks[current_metadata_idx + 1]
            if sequence_count >= next_metadata.start_index:
                current_metadata = next_metadata
                current_inverse_table = inverse_tables[current_metadata_idx + 1]
                current_metadata_idx += 1

        header = None

        if sequence_count in headers_cache:  # Try cache first
            header = headers_cache[sequence_count]
            if header and not header.startswith("@"):
                header = "@" + header

        if not header:  # Fallback to generation @seq header
            if current_metadata and current_metadata.structure_template:
                header_prefix = find_structure_prefix(
                    current_metadata.structure_template
                )
                header = (
                    f"@{header_prefix}.{sequence_count}"
                    if header_prefix
                    else f"@seq{sequence_count}"
                )
            elif sra_accession:
                header = f"@{sra_accession}.{sequence_count}"
            else:
                header = f"@seq{sequence_count}"
            logger.warning(f"Using fallback header for seq {sequence_count}: {header}")

        seq_array = np.frombuffer(seq_data, dtype=np.uint8)
        bases_array = reverse_map[seq_array]
        bases = bases_array.tobytes().decode("ascii")

        # Reconstruct quality
        if current_inverse_table is not None:
            quality_scores = reverse_scaling_to_quality(
                seq_array,
                subtract_table,
                current_inverse_table,
                max_phred=phred_alphabet_max,
            )
            quality_string = bytes(
                (quality_scores + phred_offset).astype(np.uint8)
            ).decode("ascii")

        # Write output
        output_buffer.write(header.encode("utf-8"))
        output_buffer.write(b"\n")
        output_buffer.write(bases.encode("ascii"))

        plus_line = b"\n+"
        if second_head_flag and header:
            plus_line += header[1:].encode("utf-8")
        plus_line += b"\n"
        output_buffer.write(plus_line)

        output_buffer.write(quality_string.encode("ascii"))
        output_buffer.write(b"\n")

        sequence_count += 1
        sequences_in_chunk += 1

    return (output_buffer.getvalue(), sequences_in_chunk)
