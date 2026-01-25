import argparse
import cProfile
import logging
import os
import pstats
import time
from functools import partial
from multiprocessing import Pool

import numpy as np

from toFASTR_chunk_processor import chunk_generator, process_chunk_worker
from toFASTR_header_compression import (format_metadata_header,
                                        metadata_dict_equals)
from toFASTR_quality_processing import (create_phred_quality_map,
                                        get_scaling_equation,
                                        validate_and_adjust_formula)

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

__version__ = "1.0.0"


def convert_fastq_to_fastr(
    fastq_path,
    base_map,
    output_path,
    phred_map=None,
    min_quality=0,
    quality_scaling="none",
    binary=True,
    compress_headers=False,
    sequencer_type="none",
    sra_accession=None,
    keep_bases=False,
    keep_quality=False,
    custom_formula=None,
    multiple_flowcells=False,
    remove_repeating_header=False,
    phred_alphabet_max=41,
    paired_end=False,
    paired_end_mode="same_file",
    chunk_size_mb=8,
    num_workers=4,
    adaptive_sample_size=10,
    mode=None,
    mode3_input_headers=None,
    second_head=None,
    safe_mode=False,
    verbose=False,
):
    """
    Convert FASTQ files to FASTR compressed format by encoding bases as integers and optionally compressing repetitive header metadata.
    """

    file_size = os.path.getsize(fastq_path)
    logger.info(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    logger.info(f"Using {num_workers} worker processes for parallel processing")

    # Are we extracting/reading headers for mode 3?
    extract_headers = mode == 3 and mode3_input_headers is None
    read_headers = mode == 3 and mode3_input_headers is not None
    if extract_headers:
        headers_output_path = output_path.rsplit(".", 1)[0] + "_headers.txt"
        logger.info(f"Mode 3: Extracting headers to {headers_output_path}")

    elif read_headers:
        logger.info(f"Mode 3: Reading headers from {mode3_input_headers}")

    if compress_headers and sequencer_type != "none":
        logger.info(f"Header compression enabled (sequencer type: {sequencer_type})")
        if multiple_flowcells:
            logger.info("Multiple flowcell detection enabled")
        if remove_repeating_header:
            logger.info(
                "Repeating header metadata removal enabled - only unique IDs will be stored"
            )

    if paired_end:
        logger.info(f"Paired-end mode enabled (output mode: {paired_end_mode})")

    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    logger.info(f"Processing with {chunk_size_mb}MB chunks")

    # Premake lookup table
    BYTE_LOOKUP = np.array([f"{i:03d}".encode("ascii") for i in range(256)])

    # For tracking flowcells and sequences
    flowcell_metadata_list = []
    current_flowcell_metadata = None
    structure_template = None
    delimiter = None
    flowcell_start_index = 0
    total_sequences = 0

    headers_file = None
    if extract_headers:
        headers_file = open(headers_output_path, "wb", buffering=chunk_size_bytes)
    elif read_headers:
        headers_file = open(mode3_input_headers, "rb")
        headers_file.close()

    with open(output_path, "wb", buffering=chunk_size_bytes) as outfile:
        # Write SRA accession if provided
        if sra_accession:
            outfile.write(f"#{sra_accession}\n".encode("utf-8"))

        # Create chunk generator
        chunk_gen = chunk_generator(fastq_path, chunk_size_bytes)
        if verbose:
            logger.debug("Starting chunk generation and processing...")

        first_chunk_data = next(
            chunk_gen
        )  # Process FIRST chunk to get structure/delimiter
        if verbose:
            logger.info(f"Processing first chunk to detect structure/delimiter...")
        worker_func_first = (
            partial(  # Create worker function WITHOUT structure initially
                process_chunk_worker,
                base_map=base_map,
                phred_map=phred_map,
                compress_headers=compress_headers,
                sequencer_type=sequencer_type,
                keep_bases=keep_bases,
                keep_quality=keep_quality,
                quality_scaling=quality_scaling,
                custom_formula=custom_formula,
                phred_alphabet_max=phred_alphabet_max,
                min_quality=min_quality,
                BYTE_LOOKUP=BYTE_LOOKUP,
                binary=binary,
                remove_repeating_header=remove_repeating_header,
                adaptive_structure=None,
                adaptive_delimiter=None,
                adaptive_sample_size=adaptive_sample_size,
                extract_headers=extract_headers,
                mode=mode,
                safe_mode=safe_mode,
                verbose=verbose,
            )
        )
        (
            chunk_id,
            first_processed_bytes,
            metadata,
            structure_template,
            delimiter,
            count,
            first_headers_data,
        ) = worker_func_first(first_chunk_data)

        total_sequences = count

        if metadata:
            current_flowcell_metadata = metadata
            flowcell_start_index = 0

        if sequencer_type != "none" or mode is not None:
            metadata_lines = []
            metadata_lines.append(f"#MODE={mode}\n")

            seq_type_display = sequencer_type
            if sra_accession:
                if "_sra" not in sequencer_type:
                    seq_type_display = f"{sequencer_type}_sra"

            metadata_lines.append(f"#SEQ-TYPE={seq_type_display}\n")

            # For now, write single flowcell metadata (we'll handle multiple flowcells in a second pass if needed)
            metadata_obj = (
                current_flowcell_metadata if current_flowcell_metadata else {}
            )
            write_sequencer_metadata(
                metadata_lines,
                metadata_obj,
                sequencer_type,
                structure_template,
                paired_end,
                paired_end_mode,
                quality_scaling,
                custom_formula,
                phred_alphabet_max,
                0,
                -1,
                base_map,
                mode,
                second_head,
                safe_mode,
            )
            if verbose:
                logger.info(f"First chunk processed: {count} sequences")
                if structure_template:
                    logger.info(f"Detected structure: {structure_template}")
                if delimiter:
                    logger.info(f"Detected delimiter: '{delimiter}'")

            metadata_bytes = "".join(metadata_lines).encode("utf-8")
            outfile.write(metadata_bytes)

        outfile.write(first_processed_bytes)
        if extract_headers and first_headers_data:
            headers_file.write(first_headers_data)

        # Create worker pool and stream chunks as they're processed
        with Pool(processes=num_workers) as pool:
            # Create partial function w/ fixed args
            worker_func = partial(
                process_chunk_worker,
                base_map=base_map,
                phred_map=phred_map,
                compress_headers=compress_headers,
                sequencer_type=sequencer_type,
                keep_bases=keep_bases,
                keep_quality=keep_quality,
                quality_scaling=quality_scaling,
                custom_formula=custom_formula,
                phred_alphabet_max=phred_alphabet_max,
                min_quality=min_quality,
                BYTE_LOOKUP=BYTE_LOOKUP,
                binary=binary,
                remove_repeating_header=remove_repeating_header,
                adaptive_structure=structure_template,
                adaptive_delimiter=delimiter,
                adaptive_sample_size=adaptive_sample_size,
                extract_headers=extract_headers,
                mode=mode,
                safe_mode=safe_mode,
                verbose=verbose,
            )
            chunk_counter = 1
            for (
                chunk_id,
                processed_bytes,
                metadata,
                structure,
                delimiter_result,
                count,
                headers_data,
            ) in pool.imap(worker_func, chunk_gen, chunksize=1):
                if verbose:
                    logger.debug(
                        f"Processing chunk {chunk_counter}: {count} sequences, {len(processed_bytes)} bytes written"
                    )
                    chunk_counter += 1
                outfile.write(processed_bytes)
                if extract_headers and headers_data:
                    headers_file.write(headers_data)

                # Capture structure template from first chunk
                if structure and structure_template is None:
                    structure_template = structure
                # Capture delimiter from first chunk
                if delimiter_result and delimiter is None:
                    delimiter = delimiter_result

                # Track flowcell metadata
                if metadata:
                    if multiple_flowcells:
                        if current_flowcell_metadata is None:
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                        elif not metadata_dict_equals(
                            current_flowcell_metadata, metadata
                        ):
                            flowcell_metadata_list.append(
                                (
                                    current_flowcell_metadata.copy(),
                                    flowcell_start_index,
                                    total_sequences - 1,
                                    structure_template,
                                )
                            )
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                            logger.debug(
                                f"Flowcell change detected at sequence {total_sequences}"
                            )
                    elif current_flowcell_metadata is None:
                        current_flowcell_metadata = metadata

                total_sequences += count

        # Finalize flowcell tracking
        if multiple_flowcells and current_flowcell_metadata is not None:
            flowcell_metadata_list.append(
                (
                    current_flowcell_metadata.copy(),
                    flowcell_start_index,
                    total_sequences - 1,
                    structure_template,
                )
            )
        elif current_flowcell_metadata:
            flowcell_metadata_list.append(
                (current_flowcell_metadata, 0, total_sequences - 1, structure_template)
            )

        # Print flowcell summary
        if multiple_flowcells and len(flowcell_metadata_list) > 1:
            logger.warning(
                f"\nDetected {len(flowcell_metadata_list)} flowcells - metadata written for first flowcell only!"
            )
            logger.warning(
                "Multiple flowcell support requires two-pass processing or index-based reconstruction."
            )
            for idx, (metadata, start, end, _) in enumerate(flowcell_metadata_list):
                fc_id = metadata.get("flowcell", metadata.get("movie", "unknown"))
                logger.info(
                    f"  Flowcell {idx + 1}: {fc_id} (sequences {start}-{end}, total: {end - start + 1})"
                )

    if extract_headers and headers_file:
        headers_file.close()
        logger.info(f"Headers saved to: {headers_output_path}")

    logger.info(f"Total sequences written: {total_sequences:,}")


def write_sequencer_metadata(
    metadata_lines,
    metadata,
    sequencer_type,
    structure,
    paired_end,
    paired_end_mode,
    quality_scaling,
    custom_formula,
    phred_alphabet_max,
    start_idx,
    end_idx,
    grayscale_map,
    mode,
    second_head,
    safe_mode,
):
    """
    Write metadata headers for different sequencer types.
    """
    instrument = ""
    sample_id = ""
    run_number = ""
    run_date = ""
    run_start_time = ""
    flow_cell = ""
    accession = ""
    length = ""  # Should be 'y' if length=x exists in header, blank otherwise

    if "sra" in metadata:
        accession = metadata["sra"]
    # This cleans the accession (so we dont mess up and write out something like SRR.X in the metadata, but rather SRR)
    if accession:
        for delimiter in [".", " ", "_", "/", ":"]:
            if delimiter in accession:
                accession = accession.split(delimiter)[0]
                break

    if sequencer_type in ["pacbio_hifi_sra", "pacbio_clr_sra"]:
        movie = metadata.get("movie", "")
        if movie and "_" in movie:
            parts = movie.split("_")
            instrument = parts[0]  # e.g., 'm64011'
            if len(parts) >= 2 and len(parts[1]) == 6:
                run_date = parts[1]  # e.g., '190830'
            if len(parts) >= 3 and len(parts[2]) == 6:
                run_start_time = parts[2]  # e.g., '220126'
        else:
            instrument = movie

    elif sequencer_type == "illumina_sra":
        instrument = metadata.get("instrument", "")
        flow_cell = metadata.get("flowcell", "")
        run_number = metadata.get("run_id", "")
        if structure:
            parts = structure.split()
            if len(parts) > 1:
                header_parts = parts[1].split(":")

    elif sequencer_type == "ont_sra":
        sample_id = metadata.get("sampleid", "")
        run_number = metadata.get("runid", "")

    elif sequencer_type in ["illumina", "old_illumina"]:
        instrument = metadata.get("instrument", "")
        flow_cell = metadata.get("flowcell", "")
        run_number = metadata.get("run_id", "")

    elif sequencer_type.startswith("pacbio"):
        movie = metadata.get("movie", "")
        if movie and "_" in movie:
            parts = movie.split("_")
            instrument = parts[0]
            if len(parts) >= 2 and len(parts[1]) == 6:
                run_date = parts[1]
            if len(parts) >= 3 and len(parts[2]) == 6:
                run_start_time = parts[2]
        else:
            instrument = movie

    elif sequencer_type == "ont":
        sample_id = metadata.get("sampleid", "")
        run_number = metadata.get("runid", "")

    elif sequencer_type == "sra":
        prefix = metadata.get("prefix", "")
        accession_num = metadata.get("accession", "")
        if prefix and accession_num:
            accession = f"{prefix}{accession_num}"

    # Check if length=x exists in structure (ONT format typically has this)
    # OR if has_length flag is set in metadata (PacBio SRA formats)
    if (structure and "length=" in structure.lower()) or metadata.get(
        "has_length", False
    ):
        length = "y"
    if second_head:
        second_head = "y"

    metadata_lines.append(f"#ACCESSION={accession}\n")
    metadata_lines.append(f"#INSTRUMENT={instrument}\n")
    metadata_lines.append(f"#SAMPLE-ID={sample_id}\n")
    metadata_lines.append(f"#RUN-NUMBER={run_number}\n")
    metadata_lines.append(f"#RUN-DATE(YYYYMMDD)={run_date}\n")
    metadata_lines.append(f"#RUN-START-TIME={run_start_time}\n")
    metadata_lines.append(f"#FLOW-CELL={flow_cell}\n")
    metadata_lines.append(
        f"#PHRED-ALPHABET=PHRED_{phred_alphabet_max + 1}\n"
    )  # We add one because the raw "Max" is one less than the values it can represent in total, due to the inclusion of 0
    metadata_lines.append(
        f"#GRAY_VALS={grayscale_map[[ord('N')]]},{grayscale_map[[ord('A')]]},{grayscale_map[[ord('G')]]},{grayscale_map[[ord('C')]]},{grayscale_map[[ord('T')]]}\n"
    )
    metadata_lines.append(f"#LENGTH={length}\n")
    metadata_lines.append(f"#SECOND_HEAD={second_head}\n")
    metadata_lines.append(f"#SAFE_MODE={('y' if safe_mode else '')}\n")
    metadata_lines.append(f"#PAIRED-END={'1' if paired_end else ''}\n")
    metadata_lines.append(
        f"#PAIRED-END-SAME-FILE={'1' if (paired_end and paired_end_mode == 'same_file') else ''}\n"
    )

    equation = get_scaling_equation(quality_scaling, custom_formula, phred_alphabet_max)
    metadata_lines.append(f"#QUAL_SCALE={equation}\n")

    # Don't write structure for mode 3 since headers are in separate file
    if mode == 3:
        metadata_lines.append(f"#STRUCTURE=\n")
    elif structure:
        metadata_lines.append(f"#STRUCTURE={structure}\n")
    else:
        metadata_lines.append(f"#STRUCTURE=\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert and compress FASTQ/FASTA files to FASTR format.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Positional Arguments
    parser.add_argument(
        "input_path", metavar="FILE", help="Path of .fasta or .fastq file"
    )
    parser.add_argument("output_path", metavar="FILE", help="Output file path")

    # Mode Group
    mode_group = parser.add_argument_group("OPERATION MODES")
    mode_group.add_argument(
        "--mode",
        type=int,
        metavar="INT",
        help="0: Header compression only\n"
        "1: Base conversion into numbers only\n"
        "2: Header and base conversion, written out in two lines\n"
        "3: Repeating header removal entirely, base conversion kept, written out in one line",
    )

    # Quality Group
    quality_group = parser.add_argument_group("QUALITY SCALING")
    quality_group.add_argument(
        "--qual_scale",
        type=str,
        metavar="STR",
        choices=["log", "log_reverse", "log_custom", "one_hot", "custom"],
        default="one_hot",
        help="Quality scaling method. Available options: {'log', 'log_reverse', 'log_custom', 'one_hot', 'custom'} [one_hot]",
    )
    quality_group.add_argument(
        "--extract_qual",
        type=int,
        default=1,
        metavar="INT",
        help="For FASTQ: extract quality scores (0/1) [1]",
    )
    quality_group.add_argument(
        "--phred_off",
        type=int,
        default=33,
        metavar="INT",
        help="Phred quality offset [33]",
    )
    quality_group.add_argument(
        "--min_qual",
        type=int,
        default=0,
        metavar="INT",
        help="Clamped minimum quality score threshold [0]",
    )
    quality_group.add_argument(
        "--custom_formula",
        type=str,
        default=None,
        metavar="STR",
        help="Custom formula for quality scaling (use 'x' for quality score). "
        "Example: '1 + 62 * (x - 40) / 53' or 'ln(x) * 10'",
    )

    # Paired-end args
    paired_end = parser.add_argument_group("PAIRED-END")
    paired_end.add_argument(
        "--paired",
        type=int,
        default=0,
        metavar="INT",
        help="Paired-end reads flag (0/1) [0]",
    )
    paired_end.add_argument(
        "--paired_mode",
        type=str,
        metavar="STR",
        choices=["same_file", "separate_files"],
        default="same_file",
        help="Output mode for paired-end reads. Available options: {'same_file', 'separate_files'} [same_file]",
    )

    # Sequencer Group
    seq_group = parser.add_argument_group("SEQUENCER & HEADERS")
    seq_group.add_argument(
        "--seq_type",
        type=str,
        metavar="STR",
        choices=[
            "none",
            "adaptive",
            "illumina",
            "pacbio_ccs",
            "pacbio_hifi",
            "pacbio_subread",
            "pacbio_clr",
            "ont",
            "sra",
            "old_illumina",
            "pacbio_hifi_sra",
            "pacbio_clr_sra",
            "ont_sra",
            "illumina_sra",
        ],
        default="adaptive",
        help=(
            "Sequencer type for header compression. [adaptive]\n"
            "Standard: {'illumina', 'pacbio_hifi', 'pacbio_clr', 'ont', 'sra', 'old_illumina'}\n"
            "SRA Hybrid: {'illumina_sra', 'pacbio_hifi_sra', 'pacbio_clr_sra', 'ont_sra'}"
        ),
    )
    seq_group.add_argument(
        "--compress_hdr",
        type=int,
        default=0,
        metavar="INT",
        help="Compress FASTQ headers on-the-fly (0/1) [0]",
    )
    seq_group.add_argument(
        "--sra_acc",
        type=str,
        default=None,
        metavar="STR",
        help="SRA accession number (e.g., SRR12345678) [null]",
    )
    seq_group.add_argument(
        "--multi_flow",
        type=int,
        default=0,
        metavar="INT",
        help="Enable multiple flowcell detection and tracking (0/1) [0]",
    )
    seq_group.add_argument(
        "--rm_repeat_hdr",
        type=int,
        default=0,
        metavar="INT",
        help="Remove repeating metadata from headers, store only at top (0/1) [0]",
    )
    seq_group.add_argument(
        "--adaptive_sample",
        type=int,
        default=10,
        metavar="INT",
        help="Number of headers to analyze for adaptive pattern detection [10]",
    )
    parser.add_argument(
        "--mode3_headers",
        type=str,
        default=None,
        metavar="STR",
        help="Path to headers file for mode 3 reconstruction (read mode) [null]",
    )

    # Encoding/Grayscale Group
    gray_group = parser.add_argument_group("ENCODING & GRAYSCALE")
    gray_group.add_argument(
        "--gray_N", type=int, default=0, metavar="INT", help="Grayscale value for N [0]"
    )
    gray_group.add_argument(
        "--gray_A", type=int, default=3, metavar="INT", help="Grayscale value for A [3]"
    )
    gray_group.add_argument(
        "--gray_G",
        type=int,
        default=66,
        metavar="INT",
        help="Grayscale value for G [66]",
    )
    gray_group.add_argument(
        "--gray_C",
        type=int,
        default=129,
        metavar="INT",
        help="Grayscale value for C [129]",
    )
    gray_group.add_argument(
        "--gray_T",
        type=int,
        default=192,
        metavar="INT",
        help="Grayscale value for T [192]",
    )

    # Output format
    output_group = parser.add_argument_group("OUTPUT FORMAT")
    output_group.add_argument(
        "--bin_write",
        type=int,
        default=1,
        metavar="INT",
        help="Enable binary writing of sequence integers (0/1) [1]",
    )
    output_group.add_argument(
        "--keep_bases",
        type=int,
        default=0,
        metavar="INT",
        help="Return textual bases without scaling or one-hot encoding (0/1) [0]",
    )
    output_group.add_argument(
        "--keep_qual",
        type=int,
        default=0,
        metavar="INT",
        help="Keep original quality scores in output (0/1) [0]",
    )
    output_group.add_argument(
        "--phred_alpha",
        type=str,
        default="phred94",
        metavar="STR",
        help="Phred quality (q-score) ascii character alphabet used by input (phred42, phred63, phred94) [phred94]",
    )
    output_group.add_argument(
        "--second_head",
        type=int,
        default=0,
        metavar="INT",
        help="Repeat the header on the '+' line in the FASTQ output.",
    )
    output_group.add_argument(
        "--safe_mode",
        type=int,
        default=1,
        metavar="INT",
        help="Enable safe mode for modes 1 and 2 (adds 255 marker after headers) (0/1) [1]",
    )

    # Performance Group
    perf_group = parser.add_argument_group("PERFORMANCE & PARALLELIZATION")
    perf_group.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="INT",
        help="Number of parallel workers (use 4+ for large files >5GB) [1]",
    )
    perf_group.add_argument(
        "--chunk_mb",
        type=int,
        default=8,
        metavar="INT",
        help="Chunk size in MB for parallel processing [8]",
    )
    perf_group.add_argument(
        "--profile",
        type=int,
        default=0,
        metavar="INT",
        help="Enable profiling (0/1) [0]",
    )
    perf_group.add_argument(
        "--verbose",
        type=int,
        default=0,
        metavar="INT",
        help="Enable verbose logging (0/1) [0]",
    )

    args = parser.parse_args()

    if args.verbose == 1:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("toFASTR_chunk_processor").setLevel(logging.DEBUG)
        logging.getLogger("toFASTR_fastr_parser").setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if args.qual_scale == "custom" and args.custom_formula is None:
        logger.error("--qual_scale custom requires --custom_formula argument")
        logger.info("Example usage:")
        logger.info("  --qual_scale custom --custom_formula '1 + 62 * (x - 40) / 53'")
        exit(1)

    if args.multi_flow == 1 and args.compress_hdr == 0:
        logger.warning("--multi_flow requires --compress_hdr 1 to function")
        logger.info("Enabling header compression automatically...")
        args.compress_hdr = 1

    if args.rm_repeat_hdr == 1 and args.compress_hdr == 0:
        logger.warning("WARNING: --rm_repeat_hdr requires --compress_hdr 1 to function")
        logger.info("Enabling header compression automatically...")
        args.compress_hdr = 1

    # Mode shortcuts
    if args.mode == 0:
        args.compress_hdr = 1
        args.bin_write = 0
        args.keep_bases = 1
        args.keep_qual = 1
        args.safe_mode = 0  # Safe mode not needed for mode 0 due to a lack of binary (Just compressed headers and not-processed base/qualities)
    elif args.mode == 1:
        args.compress_hdr = 0
        args.bin_write = 1
    elif args.mode == 2:
        args.compress_hdr = 1
        args.bin_write = 1
    elif args.mode == 3:
        args.compress_hdr = 1
        args.rm_repeat_hdr = 1
        args.bin_write = 1
        args.safe_mode = 0  # Safe mode not needed for mode 3 (as it is already safe)

    # Phred alphabet configuration
    if args.phred_alpha == "phred42":
        phred_alphabet_max = 41
    elif args.phred_alpha == "phred63":
        phred_alphabet_max = 62
    elif args.phred_alpha == "phred94":
        phred_alphabet_max = 93

    # Validate custom formula if provided
    if args.custom_formula:
        args.custom_formula = validate_and_adjust_formula(
            args.custom_formula, phred_alphabet_max
        )

    # Create base map
    numpy_base_map = np.zeros(128, dtype=np.uint8)
    numpy_base_map[ord("N")] = args.gray_N
    numpy_base_map[ord("A")] = args.gray_A
    numpy_base_map[ord("G")] = args.gray_G
    numpy_base_map[ord("C")] = args.gray_C
    numpy_base_map[ord("T")] = args.gray_T

    max_base_value = max(
        args.gray_N, args.gray_A, args.gray_G, args.gray_C, args.gray_T
    )
    if (
        max_base_value + 62 >= 255
    ):  # 62 is max quality scaling range, so we make sure nothing meets/exceeds our reserved val (255)
        logger.error(
            f"Base encoding values too high. Maximum base value ({max_base_value}) + quality range (62) must be < 255"
        )
        logger.error(
            f"Current values: N={args.gray_N}, A={args.gray_A}, G={args.gray_G}, C={args.gray_C}, T={args.gray_T}"
        )
        exit(1)

    start_time = time.perf_counter()

    # Initialize profiler if requested
    profiler = None
    if args.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        print("Profiling enabled...")

    # Create phred map if needed
    phred_map = (
        create_phred_quality_map(args.phred_off, phred_alphabet_max)
        if args.extract_qual
        else None
    )

    logger.info(f"Converting sequences from {args.input_path}...")

    convert_fastq_to_fastr(
        args.input_path,
        numpy_base_map,
        args.output_path,
        phred_map,
        args.min_qual,
        args.qual_scale,
        args.bin_write,
        compress_headers=(args.compress_hdr == 1),
        sequencer_type=args.seq_type,
        sra_accession=args.sra_acc,
        keep_bases=args.keep_bases,
        keep_quality=args.keep_qual,
        custom_formula=args.custom_formula,
        multiple_flowcells=(args.multi_flow == 1),
        remove_repeating_header=(args.rm_repeat_hdr == 1),
        phred_alphabet_max=phred_alphabet_max,
        paired_end=args.paired,
        paired_end_mode=args.paired_mode,
        num_workers=args.workers,
        chunk_size_mb=args.chunk_mb,
        adaptive_sample_size=args.adaptive_sample,
        mode=args.mode,
        mode3_input_headers=args.mode3_headers,
        second_head=args.second_head,
        safe_mode=(args.safe_mode == 1),
        verbose=(args.verbose == 1),
    )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    logger.info(f"Conversion completed in {elapsed_time:.4f} seconds")
    logger.info(f"Output saved to {args.output_path}")

    if profiler is not None:
        profiler.disable()
        print("\n" + "=" * 80)
        print("Profiling Results:")
        print("=" * 80)
        stats = pstats.Stats(profiler)
        stats.sort_stats("cumulative")
        stats.print_stats(30)
        stats.sort_stats("tottime")
        print("\n" + "=" * 80)
        print("By total time:")
        print("=" * 80)
        stats.print_stats(30)


if __name__ == "__main__":
    main()
