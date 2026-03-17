import io
import logging
import traceback
from typing import Tuple

import numpy as np

from toFASTR_fastq_parser import parse_fastq_records_from_buffer
from toFASTR_fastr_writer import process_and_write_records

logger = logging.getLogger(__name__)


def process_chunk_worker(
    chunk_data: Tuple[int, bytes, int],
    base_map: np.ndarray,
    phred_map: np.ndarray,
    compress_headers: bool,
    sequencer_type: str,
    keep_bases: bool,
    keep_quality: bool,
    quality_scaling: str,
    custom_formula: str,
    phred_alphabet_max: int,
    min_quality: int,
    BYTE_LOOKUP: np.ndarray,
    binary: bool,
    remove_repeating_header: bool,
    adaptive_structure: str,
    adaptive_delimiter: str,
    adaptive_sample_size: int,
    extract_headers: bool = False,
    mode: int = None,
    safe_mode: bool = False,
    verbose: bool = False,
):
    """
    Parse FASTQ records from a chunk and write to FASTR format.

    Returns:
        tuple: (chunk_id, processed_bytes, metadata, structure, delimiter, count, headers_data)
    """
    try:
        chunk_id, buffer, start_index = chunk_data
        if verbose:
            logger.debug(f"Worker processing chunk {chunk_id} with {len(buffer)} bytes")

        # Parse the records from buffer
        records, _, metadata, structure, delimiter, count = (
            parse_fastq_records_from_buffer(
                buffer,
                start_index,
                base_map,
                phred_map,
                compress_headers,
                sequencer_type,
                keep_bases,
                keep_quality,
                adaptive_structure=adaptive_structure,
                adaptive_delimiter=adaptive_delimiter,
                adaptive_sample_size=adaptive_sample_size,
            )
        )

        if verbose:
            logger.debug(f"Worker parsed {count} records from chunk {chunk_id}")

        # Process records and write to in-memory buffer
        output_buffer = io.BytesIO()
        headers_buffer = io.BytesIO() if extract_headers else None

        if records:
            process_and_write_records(
                records,
                output_buffer,
                base_map,
                quality_scaling,
                custom_formula,
                phred_alphabet_max,
                min_quality,
                keep_bases,
                binary,
                keep_quality,
                remove_repeating_header,
                compress_headers,
                BYTE_LOOKUP,
                headers_buffer=headers_buffer,
                mode=mode,
                safe_mode=safe_mode,
            )

        #  print(f"Worker completed chunk {chunk_id}")
        headers_data = headers_buffer.getvalue() if headers_buffer else None
        return (
            chunk_id,
            output_buffer.getvalue(),
            metadata,
            structure,
            delimiter,
            count,
            headers_data,
        )

    except Exception as e:
        logger.error(f"Error in worker processing chunk {chunk_id}: {e}", exc_info=True)
        raise


def find_last_fastq_record_boundary(buffer: bytes) -> int:
    """
    Find the position of the last complete FASTQ record boundary in buffer.
    Returns the position after the last complete record, or -1 if no complete record found.

    A complete FASTQ record has 4 lines:
    1. @ header
    2. sequence
    3. + separator
    4. quality (same length as sequence)
    """
    lines = buffer.split(b"\n")

    # Work backwards to find the last complete record
    i = len(lines) - 1

    # Skip empty trailing line if present
    if i >= 0 and lines[i] == b"":
        i -= 1

    # We need at least 4 lines for a complete record
    while i >= 3:
        # Check if this could be the end of a quality line (line 4 of record)
        qual_line = lines[i]
        plus_line = lines[i - 1] if i >= 1 else b""
        seq_line = lines[i - 2] if i >= 2 else b""
        header_line = lines[i - 3] if i >= 3 else b""

        # Valid record structure:
        # - Header starts with @
        # - Plus line starts with +
        # - Quality length matches sequence length
        if (
            header_line.startswith(b"@")
            and plus_line.startswith(b"+")
            and len(qual_line) == len(seq_line)
        ):
            # Found a complete record! Return position after this record
            return sum(len(line) + 1 for line in lines[: i + 1])

        i -= 1

    return -1


def chunk_generator(fastq_path: str, chunk_size_bytes: int):
    """
    Generator that yields file chunks ending on complete record boundaries.

    Args:
        fastq_path: Path to FASTQ input file
        chunk_size_bytes: Size of chunks to read

    Yields:
        tuple: (chunk_id, buffer_data, start_record_index)
    """
    buffer = b""
    chunk_id = 0
    start_index = 0
    file_position = 0

    with open(fastq_path, "rb", buffering=chunk_size_bytes) as infile:
        while True:
            chunk = infile.read(chunk_size_bytes)
            if not chunk and not buffer:
                break

            buffer += chunk
            file_position += len(chunk)

            # Find last complete record boundary
            if not chunk:
                # End of file - process all remaining buffer
                process_buffer = buffer
                buffer = b""
            else:
                # Find last complete FASTQ record
                boundary_pos = find_last_fastq_record_boundary(buffer)
                if boundary_pos == -1:
                    # No complete record found, keep accumulating
                    continue
                else:
                    process_buffer = buffer[:boundary_pos]
                    buffer = buffer[boundary_pos:]

            if process_buffer:
                # The actual record count will be determined by parse_fastq_records_from_buffer
                # For now, we provide an estimate for the start_index
                # This estimate doesn't need to be perfect since records have absolute headers
                lines_in_chunk = process_buffer.count(b"\n")
                estimated_records = lines_in_chunk // 4  # FASTQ has 4 lines per record

                yield (chunk_id, process_buffer, start_index)

                logger.debug(
                    f"Yielding chunk {chunk_id}: bytes 0-{file_position} ({len(process_buffer)} bytes, ~{estimated_records} seqs)"
                )
                start_index += estimated_records
                chunk_id += 1
