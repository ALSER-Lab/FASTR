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
    Worker function that parses FASTQ records and applies transformations.

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
            last_at = buffer.rfind(b"\n@")
            if last_at == -1 or not chunk:
                process_buffer = buffer
                buffer = b""
            else:
                process_buffer = buffer[: last_at + 1]
                buffer = buffer[last_at + 1 :]

            if process_buffer:
                yield (chunk_id, process_buffer, start_index)

                # Estimate number of records for next chunk's start index
                estimated_records = process_buffer.count(b"\n@")
                logger.debug(
                    f"Yielding chunk {chunk_id}: bytes 0-{file_position} ({len(process_buffer)} bytes, ~{estimated_records} seqs)"
                )
                start_index += estimated_records
                chunk_id += 1

