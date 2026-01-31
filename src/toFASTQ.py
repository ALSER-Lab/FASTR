import argparse
import cProfile
import logging
import os
import pstats
import time
from functools import partial
from io import StringIO
from multiprocessing import Pool

import numpy as np

from toFASTQ_base_mapping import create_base_map, reverse_base_map
from toFASTQ_header_indexing import build_header_index
from toFASTQ_metadata_parser import parse_metadata_header
from toFASTQ_quality_reconstruction import (build_formula_func,
                                            build_inverse_quality_table)
from toFASTQ_reconstruction_worker import process_chunk_worker_reconstruction

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def reconstruct_fastq(input_path: str, output_path: str, **kwargs):
    """
    Main function to reconstruct FASTQ file from FASTR.
    Coordinates parallel reconstruction w/ worker processes.
    """
    chunk_size_mb = kwargs.get("chunk_size_mb", 8)
    num_workers = kwargs.get("threads", 1)
    mode3_headers_file = kwargs.get("mode3_headers_file", None)
    verbose = kwargs.get("verbose", False)

    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    file_size = os.path.getsize(input_path)

    with open(input_path, "rb") as f:
        header_data = f.read(min(1024 * 1024, file_size))

    (
        metadata_blocks,
        data_start_byte,
        sra_acc,
        phred_from_metadata,
        detected_mode,
        length_flag,
        second_head_flag,
        safe_mode_flag,
        grayscale_vals,
    ) = parse_metadata_header(header_data)
    if grayscale_vals is not None and len(grayscale_vals) == 5:
        gray_N, gray_A, gray_G, gray_C, gray_T = grayscale_vals
        logger.info(
            f"Using grayscale values from metadata: N={gray_N}, A={gray_A}, G={gray_G}, C={gray_C}, T={gray_T}"
        )
    else:
        gray_N, gray_A, gray_G, gray_C, gray_T = (
            0,
            3,
            66,
            129,
            192,
        )  # Default vals
        logger.info("Using default grayscale values")

    mode = detected_mode
    phred_alphabet_max = phred_from_metadata
    phred_offset = 33

    subtract_table = create_base_map(gray_N, gray_A, gray_G, gray_C, gray_T)
    reverse_map = reverse_base_map(gray_N, gray_A, gray_G, gray_C, gray_T)

    inverse_tables = []
    if metadata_blocks:
        for mb in metadata_blocks:
            formula_func = build_formula_func(mb.scaling_equation)

            q_possible = np.arange(phred_alphabet_max + 1, dtype=np.float32)
            scaled_int_for_q = formula_func(q_possible).astype(np.int32)

            inverse_table = build_inverse_quality_table(
                scaled_int_for_q, q_possible.astype(np.int32), 63
            )
        inverse_tables.append(inverse_table)
    headers_mmap_info = None
    header_index_mmap = None
    mmap_cleanup_path = None

    if mode == 3 and mode3_headers_file:
        header_index_mmap, mmap_cleanup_path = build_header_index(mode3_headers_file)
        headers_mmap_info = (
            mmap_cleanup_path,
            header_index_mmap.shape,
            header_index_mmap.dtype,
        )

    def chunk_generator():
        chunk_id = 0
        start_seq_idx = 0

        if mode == 0:
            logger.info("Mode 0: Building record index...")
            record_positions = [data_start_byte]

            with open(input_path, "rb") as f:
                f.seek(data_start_byte)
                SCAN_BLOCK_SIZE = 50 * 1024 * 1024
                buffer = b""
                file_pos = data_start_byte

                while True:
                    chunk = f.read(SCAN_BLOCK_SIZE)
                    if not chunk:
                        break

                    buffer += chunk
                    cursor = 0
                    # We run off the assumption each fastq record is exactly 4 lines, so we'll count...
                    while cursor < len(buffer):

                        if buffer[cursor : cursor + 1] != b"@":
                            cursor += 1
                            continue
                        line1_end = buffer.find(b"\n", cursor)
                        if line1_end == -1:
                            break

                        line2_start = line1_end + 1
                        line2_end = buffer.find(b"\n", line2_start)
                        if line2_end == -1:
                            break

                        bases_len = line2_end - line2_start
                        line3_start = line2_end + 1
                        line3_end = buffer.find(b"\n", line3_start)
                        if line3_end == -1:
                            break

                        if buffer[line3_start : line3_start + 1] != b"+":
                            cursor += 1
                            continue

                        line4_start = line3_end + 1
                        line4_end = line4_start + bases_len

                        if line4_end > len(buffer):
                            break

                        cursor = line4_end
                        if (
                            cursor < len(buffer)
                            and buffer[cursor : cursor + 1] == b"\n"
                        ):
                            cursor += 1
                        if cursor < len(buffer):
                            record_positions.append(file_pos + cursor)

                    buffer = buffer[cursor:]
                    file_pos += cursor

            logger.info(f"Mode 0: Found {len(record_positions)} record boundaries")

            target_bytes_per_chunk = chunk_size_bytes
            current_chunk_start = 0
            current_chunk_byte_count = 0

            for i in range(len(record_positions)):
                if i == 0:
                    continue

                record_start = record_positions[i - 1]
                record_end = (
                    record_positions[i] if i < len(record_positions) - 1 else file_size
                )
                record_bytes = record_end - record_start

                if (
                    current_chunk_byte_count + record_bytes > target_bytes_per_chunk
                    and current_chunk_byte_count > 0
                ):

                    abs_start = record_positions[current_chunk_start]
                    abs_end = record_positions[i - 1]
                    seqs_in_chunk = (i - 1) - current_chunk_start

                    if verbose:
                        logger.info(
                            f"Yielding chunk {chunk_id}: {abs_start}-{abs_end} "
                            + f"({abs_end-abs_start:,} bytes, {seqs_in_chunk} seqs)"
                        )

                    yield (
                        chunk_id,
                        abs_start,
                        abs_end,
                        start_seq_idx,
                    )  # If adding this record would exceed chunk size, yield current chunk

                    start_seq_idx += seqs_in_chunk
                    chunk_id += 1
                    current_chunk_start = i - 1
                    current_chunk_byte_count = record_bytes
                else:
                    current_chunk_byte_count += record_bytes

            if current_chunk_start < len(record_positions):
                abs_start = record_positions[current_chunk_start]
                abs_end = file_size
                seqs_in_chunk = len(record_positions) - current_chunk_start

                if verbose:
                    logger.info(
                        f"Yielding final chunk {chunk_id}: {abs_start}-{abs_end} "
                        + f"({abs_end-abs_start:,} bytes, {seqs_in_chunk} seqs)"
                    )

                yield (
                    chunk_id,
                    abs_start,
                    abs_end,
                    start_seq_idx,
                )  # Yield the final chunk

            return

        # For modes 1, 2, 3: Stream through file finding chunk boundaries
        buffer_start = data_start_byte
        buffer_len = file_size - data_start_byte
        pos = 0

        with open(input_path, "rb") as f:
            while pos < buffer_len:
                chunk_end = min(pos + chunk_size_bytes, buffer_len)

                # Read search region for boundary detection
                search_start = max(pos, chunk_end - 1000)
                search_size = chunk_end - search_start + 1000  # Extra buffer for safety

                f.seek(buffer_start + search_start)
                search_region = f.read(min(search_size, buffer_len - search_start))

                # Find boundary based on mode
                if mode == 3:
                    # Mode 3 splits on 255 markers
                    last_marker = search_region.rfind(
                        b"\n\xff", 0, chunk_end - search_start
                    )

                    if last_marker != -1:
                        chunk_end = search_start + last_marker + 1

                elif mode == 2 and safe_mode_flag:
                    # Safe mode splits on validated @header\xff pattern
                    candidate_pos = min(
                        len(search_region) - 1, chunk_end - search_start
                    )
                    found_valid = False

                    while candidate_pos > 0:
                        at_pos = search_region.rfind(b"@", 0, candidate_pos)
                        if at_pos == -1:
                            break

                        # Check if this @ has \xff after it (within reasonable distance)
                        xff_pos = search_region.find(
                            b"\xff", at_pos, min(at_pos + 200, len(search_region))
                        )
                        intervening_at = search_region.find(
                            b"\n@", at_pos + 1, min(at_pos + 200, len(search_region))
                        )

                        # Valid only if we found \xff AND no \n@ comes before it
                        if xff_pos != -1 and (
                            intervening_at == -1 or xff_pos < intervening_at
                        ):
                            chunk_end = search_start + at_pos
                            found_valid = True
                            break

                        candidate_pos = at_pos

                    if not found_valid:
                        chunk_end = min(pos + chunk_size_bytes, buffer_len)

                else:
                    # Mode 1 and 2 (unsafe) split on \n@ headers
                    last_header = search_region.rfind(
                        b"\n@", 0, chunk_end - search_start
                    )

                    if last_header != -1:
                        chunk_end = search_start + last_header + 1

                if chunk_end <= pos:
                    chunk_end = min(pos + chunk_size_bytes, buffer_len)

                abs_start = buffer_start + pos
                abs_end = buffer_start + chunk_end

                if abs_end > abs_start:
                    # Estimate sequences for this chunk
                    f.seek(abs_start)
                    sample = f.read(
                        min(abs_end - abs_start, 100000)
                    )  # Sample first 100KB

                    if mode == 3:
                        estimated_seqs = sample.count(b"\xff")
                        if len(sample) < abs_end - abs_start:
                            # Extrapolate
                            estimated_seqs = int(
                                estimated_seqs * (abs_end - abs_start) / len(sample)
                            )
                    else:
                        estimated_seqs = sample.count(b"\n@")
                        if sample.startswith(b"@"):
                            estimated_seqs += 1
                        if len(sample) < abs_end - abs_start:
                            estimated_seqs = int(
                                estimated_seqs * (abs_end - abs_start) / len(sample)
                            )

                    if verbose:
                        logger.info(
                            f"Yielding chunk {chunk_id}: {abs_start}-{abs_end} ({abs_end-abs_start} bytes, ~{estimated_seqs} seqs)"
                        )

                    yield (chunk_id, abs_start, abs_end, start_seq_idx)

                    start_seq_idx += estimated_seqs
                    chunk_id += 1

                pos = chunk_end
                if pos >= buffer_len:
                    break

    logger.info("Reconstructing FASTQ...")

    chunks_list = []
    current_seq_idx = 0

    for chunk_data in chunk_generator():
        chunk_id, abs_start, abs_end, _ = chunk_data

        with open(input_path, "rb") as f:
            f.seek(abs_start)
            chunk_size = abs_end - abs_start
            chunk_bytes = f.read(chunk_size)
            actual_count = chunk_bytes.count(b"\xff")

        chunks_list.append((chunk_id, abs_start, abs_end, current_seq_idx))
        current_seq_idx += actual_count

    total_sequences = 0

    try:
        if mode == 3 and mode3_headers_file is None:
            logger.warning(
                "WARNING: No header file specified for mode 3, using fallback @seq structure"
            )

        with open(output_path, "wb", buffering=chunk_size_bytes) as outfile:
            with Pool(processes=num_workers) as pool:
                worker_func = partial(
                    process_chunk_worker_reconstruction,
                    mmap_path=input_path,
                    reverse_map=reverse_map,
                    subtract_table=subtract_table,
                    metadata_blocks=metadata_blocks,
                    inverse_tables=inverse_tables,
                    phred_alphabet_max=phred_alphabet_max,
                    phred_offset=phred_offset,
                    sra_acc=sra_acc,
                    mode=mode,
                    headers_file_path=mode3_headers_file,
                    headers_mmap_info=headers_mmap_info,
                    data_start_byte=data_start_byte,
                    length_flag=length_flag,
                    second_head_flag=second_head_flag,
                    safe_mode_flag=safe_mode_flag,
                )

                for chunk_id, processed_bytes, count in pool.imap(
                    worker_func, chunks_list, chunksize=1
                ):
                    outfile.write(processed_bytes)
                    total_sequences += count

    finally:
        if header_index_mmap is not None:
            del header_index_mmap  # Close mmap
        if mmap_cleanup_path and os.path.exists(mmap_cleanup_path):
            try:
                os.unlink(mmap_cleanup_path)
                logger.info("Cleaned up temporary index file")
            except:
                pass

    logger.info(f"\nTotal sequences reconstructed: {total_sequences:,}")
    logger.info(f"Output saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct FASTQ files from FASTR.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Positional Arguments
    parser.add_argument("input_path", metavar="FILE", help="Path to FASTR file")
    parser.add_argument("output_path", metavar="FILE", help="Output FASTQ file path")

    # Mode Group
    mode_group = parser.add_argument_group("RECONSTRUCTION MODE")
    mode_group.add_argument(
        "--headers_file",
        type=str,
        metavar="FILE",
        default=None,
        help="Path to headers file for mode 3 reconstruction [null]",
    )

    # Performance Group
    perf_group = parser.add_argument_group("PERFORMANCE & PARALLELIZATION")
    perf_group.add_argument(
        "--chunk_size_mb",
        type=int,
        metavar="INT",
        default=8,
        help="Chunk size in MB for parallel processing [8]",
    )
    perf_group.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        default=1,
        help="Number of parallel threads [1]",
    )
    perf_group.add_argument(
        "--verbose",
        type=int,
        metavar="INT",
        default=0,
        choices=[0, 1],
        help="Enable verbose logging (0/1) [0]",
    )
    perf_group.add_argument(
        "--profile",
        type=int,
        metavar="INT",
        default=0,
        choices=[0, 1],
        help="Enable cProfile profiling (0/1) [0]",
    )

    args = parser.parse_args()

    if args.verbose == 1:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    start_time = time.perf_counter()

    if args.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        logger.info("Profiling enabled...")

    reconstruct_fastq(
        args.input_path,
        args.output_path,
        chunk_size_mb=args.chunk_size_mb,
        num_workers=args.threads,
        mode3_headers_file=args.headers_file,
        verbose=(args.verbose == 1),
    )
    if args.profile == 1:
        profiler.disable()
        s = StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
        ps.print_stats(20)
        print("\n" + "=" * 80)
        print("Profiling Results:")
        print("=" * 80)
        print(s.getvalue())

    end_time = time.perf_counter()
    logger.info(f"\nReconstruction completed in {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    main()
