import os
import numpy as np
import time
import argparse
import cProfile
import pstats
from multiprocessing import Pool
from functools import partial

from header_compression import format_metadata_header, metadata_dict_equals
from quality_processing import(create_phred_quality_map, 
                               get_scaling_equation,
                               validate_and_adjust_formula)
from chunk_processor import process_chunk_worker, chunk_generator


def export_scalars_to_txt(fastq_path, base_map, output_path, phred_map=None, min_quality=0,
                          quality_scaling='none', binary=True, compress_headers=False,
                          sequencer_type='none', sra_accession=None, keep_bases=False, keep_quality=False,
                          custom_formula=None, multiple_flowcells=False,
                          remove_repeating_header=False, phred_alphabet_max=41, paired_end=False,
                          paired_end_mode='same_file', chunk_size_mb=32, num_workers=4):
    """
    Main processing function using streaming architecture to handle files of any size.
    Uses streaming with pool.imap to maintain consistent memory usage.
    """
    
    file_size = os.path.getsize(fastq_path)
    print(f"File size: {file_size:,} bytes ({file_size / (1024**3):.2f} GB)")
    print(f"Using {num_workers} worker processes for parallel processing")
    
    if compress_headers and sequencer_type != 'none':
        print(f"Header compression enabled (sequencer type: {sequencer_type})")
        if multiple_flowcells:
            print("Multiple flowcell detection enabled")
        if remove_repeating_header:
            print("Repeating header metadata removal enabled - only unique IDs will be stored")
    
    if paired_end:
        print(f"Paired-end mode enabled (output mode: {paired_end_mode})")
    
    chunk_size_bytes = chunk_size_mb * 1024 * 1024
    print(f"Processing with {chunk_size_mb}MB chunks (parallel streaming mode)")
    
    # Premake lookup table
    BYTE_LOOKUP = np.array([f'{i:03d}'.encode('ascii') for i in range(256)])
    
    # For tracking flowcells and sequences
    flowcell_metadata_list = []
    current_flowcell_metadata = None
    flowcell_start_index = 0
    total_sequences = 0
    
    print("Reading and processing chunks in parallel...")
    
    with open(output_path, 'wb', buffering=chunk_size_bytes) as outfile:
        # Write SRA accession if provided
        if sra_accession:
            outfile.write(f"#{sra_accession}\n".encode('utf-8'))
        
        # Reserve space for metadata headers (we'll come back to write these)
        metadata_position = outfile.tell()
        MAX_METADATA_LINES = 20
        placeholder_size = MAX_METADATA_LINES * 200  # 200 chars per line
        if sequencer_type != 'none':
            outfile.write(b' ' * placeholder_size)  # Write spaces as placeholder
            outfile.write(b'\n')
        
        # Create worker pool and process chunks as they come in
        with Pool(processes=num_workers) as pool:
            # Create partial function with fixed args
            worker_func = partial(
                process_chunk_worker,
                base_map=base_map,
                phred_map=phred_map,
                compress_headers=compress_headers,
                sequencer_type=sequencer_type,
                paired_end=paired_end,
                keep_bases=keep_bases,
                keep_quality=keep_quality,
                quality_scaling=quality_scaling,
                custom_formula=custom_formula,
                phred_alphabet_max=phred_alphabet_max,
                min_quality=min_quality,
                BYTE_LOOKUP=BYTE_LOOKUP,
                binary=binary,
                remove_repeating_header=remove_repeating_header
            )
            
            # Use imap to process chunks as they're read (streaming)
            # chunksize=1 ensures order is preserved
            for chunk_id, processed_bytes, metadata, count in pool.imap(
                worker_func, 
                chunk_generator(fastq_path, chunk_size_bytes), 
                chunksize=1
            ):
                # Write immediately as each chunk completes
                outfile.write(processed_bytes)
                
                # Track flowcell metadata
                if metadata:
                    if multiple_flowcells:
                        if current_flowcell_metadata is None:
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                        elif not metadata_dict_equals(current_flowcell_metadata, metadata):
                            flowcell_metadata_list.append((
                                current_flowcell_metadata.copy(),
                                flowcell_start_index,
                                total_sequences - 1
                            ))
                            current_flowcell_metadata = metadata
                            flowcell_start_index = total_sequences
                            print(f"Flowcell change detected at sequence {total_sequences}")
                    elif current_flowcell_metadata is None:
                        current_flowcell_metadata = metadata
                
                total_sequences += count
                
                if total_sequences % 100000 == 0:
                    print(f"Processed {total_sequences:,} sequences...")
        
        # Finalize flowcell tracking
        if multiple_flowcells and current_flowcell_metadata is not None:
            flowcell_metadata_list.append((
                current_flowcell_metadata.copy(),
                flowcell_start_index,
                total_sequences - 1
            ))
        elif current_flowcell_metadata:
            flowcell_metadata_list.append((current_flowcell_metadata, 0, total_sequences - 1))
        
        # Print flowcell summary
        if multiple_flowcells and len(flowcell_metadata_list) > 1:
            print(f"\nDetected {len(flowcell_metadata_list)} flowcells:")
            for idx, (metadata, start, end) in enumerate(flowcell_metadata_list):
                fc_id = metadata.get('flowcell', metadata.get('movie', 'unknown'))
                print(f"  Flowcell {idx + 1}: {fc_id} (sequences {start}-{end}, total: {end - start + 1})")
        
        # Go back and write metadata headers at the beginning
        if sequencer_type != 'none':
            end_position = outfile.tell()
            outfile.seek(metadata_position)
            
            metadata_lines = []
            
            if multiple_flowcells and len(flowcell_metadata_list) > 0:
                # Build metadata lines for multiple flowcells
                for metadata, start_idx, end_idx in flowcell_metadata_list:
                    metadata_line = format_metadata_header(metadata, sequencer_type)
                    metadata_lines.append(f"#COMMON:{metadata_line}\n")
                    metadata_lines.append(f"#SEQUENCER:{sequencer_type}\n")
                    equation = get_scaling_equation(quality_scaling, custom_formula, phred_alphabet_max)
                    metadata_lines.append(f"#{equation}\n")
                    metadata_lines.append(f"#RANGE:{start_idx}-{end_idx}\n")
            else:
                # Single flowcell
                if current_flowcell_metadata:
                    metadata_line = format_metadata_header(current_flowcell_metadata, sequencer_type)
                    metadata_lines.append(f"#COMMON:{metadata_line}\n")
                metadata_lines.append(f"#SEQUENCER:{sequencer_type}\n")
                equation = get_scaling_equation(quality_scaling, custom_formula, phred_alphabet_max)
                metadata_lines.append(f"#{equation}\n")
            
            # Write metadata and pad remaining space
            metadata_bytes = ''.join(metadata_lines).encode('utf-8')
            outfile.write(metadata_bytes)
            remaining = placeholder_size - len(metadata_bytes)
            if remaining > 0:
                outfile.write(b' ' * remaining)
            outfile.write(b'\n')
            
            outfile.seek(end_position)  # Return to end
    
    print(f"Total sequences written: {total_sequences:,}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert and compress FASTQ/FASTA files to scalar format.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Positional Arguments
    parser.add_argument("input_path", 
                        metavar="FILE",
                        help="Path of .fasta or .fastq file")
    parser.add_argument("output_path", 
                        metavar="FILE",
                        help="Output file path")

    # Mode Group
    mode_group = parser.add_argument_group("OPERATION MODES")
    mode_group.add_argument("--mode", type=int, metavar="INT",
                            help="0: Header compression only\n"
                                 "1: Base conversion into numbers only\n"
                                 "2: Header and base conversion, written out in two lines\n"
                                 "3: Repeating header removal entirely, base conversion kept, written out in one line")

    # Quality Group
    quality_group = parser.add_argument_group("QUALITY SCALING")
    quality_group.add_argument("--qual_scale", type=str, metavar="STR",
                               choices=['log', 'log_reverse', 'log_custom', 'one_hot', 'custom'], 
                               default='one_hot',
                               help="Quality scaling method. Available options: {'log', 'log_reverse', 'log_custom', 'one_hot', 'custom'} [one_hot]")
    quality_group.add_argument("--extract_qual", type=int, default=1, metavar="INT",
                               help="For FASTQ: extract quality scores (0/1) [1]")
    quality_group.add_argument("--phred_off", type=int, default=33, metavar="INT",
                               help="Phred quality offset [33]")
    quality_group.add_argument("--min_qual", type=int, default=0, metavar="INT",
                               help="Clamped minimum quality score threshold [0]")
    quality_group.add_argument("--custom_formula", type=str, default=None, metavar="STR",
                               help="Custom formula for quality scaling (use 'x' for quality score). "
                                    "Example: '1 + 62 * (x - 40) / 53' or 'ln(x) * 10'")

    # Paired-end args
    paired_end = parser.add_argument_group("PAIRED-END")
    paired_end.add_argument("--paired", type=int, default=0, metavar="INT",
                            help="Paired-end reads flag (0/1) [0]")
    paired_end.add_argument("--paired_mode", type=str, metavar="STR",
                            choices=['same_file', 'separate_files'], 
                            default='same_file',
                            help="Output mode for paired-end reads. Available options: {'same_file', 'separate_files'} [same_file]")

    # Sequencer Group
    seq_group = parser.add_argument_group("SEQUENCER & HEADERS")
    seq_group.add_argument("--seq_type", type=str, metavar="STR",
                           choices=['none', 'illumina', 'pacbio_ccs', 'pacbio_hifi', 'pacbio_subread', 'pacbio_clr', 'ont', 'srr', 'old_illumina'],
                           default='none',
                           help="Sequencer type for header compression. Available options: {'none', 'illumina', 'pacbio_ccs', 'pacbio_hifi', 'pacbio_subread', 'pacbio_clr', 'ont', 'srr', 'old_illumina'} [none]")
    seq_group.add_argument("--compress_hdr", type=int, default=0, metavar="INT",
                           help="Compress FASTQ headers on-the-fly (0/1) [0]")
    seq_group.add_argument("--sra_acc", type=str, default=None, metavar="STR",
                           help="SRA accession number (e.g., SRR12345678) [null]")
    seq_group.add_argument("--multi_flow", type=int, default=0, metavar="INT",
                           help="Enable multiple flowcell detection and tracking (0/1) [0]")
    seq_group.add_argument("--rm_repeat_hdr", type=int, default=0, metavar="INT",
                           help="Remove repeating metadata from headers, store only at top (0/1) [0]")

    # Encoding/Grayscale Group
    gray_group = parser.add_argument_group("ENCODING & GRAYSCALE")
    gray_group.add_argument("--gray_N", type=int, default=1, metavar="INT",
                            help="Grayscale value for N [1]")
    gray_group.add_argument("--gray_A", type=int, default=63, metavar="INT",
                            help="Grayscale value for A [63]")
    gray_group.add_argument("--gray_C", type=int, default=191, metavar="INT",
                            help="Grayscale value for C [191]")
    gray_group.add_argument("--gray_G", type=int, default=255, metavar="INT",
                            help="Grayscale value for G [255]")
    gray_group.add_argument("--gray_T", type=int, default=127, metavar="INT",
                            help="Grayscale value for T [127]")

    # Output format
    output_group = parser.add_argument_group("OUTPUT FORMAT")
    output_group.add_argument("--bin_write", type=int, default=1, metavar="INT",
                              help="Enable binary writing of sequence integers (0/1) [1]")
    output_group.add_argument("--keep_bases", type=int, default=0, metavar="INT",
                              help="Return textual bases without scaling or one-hot encoding (0/1) [0]")
    output_group.add_argument("--keep_qual", type=int, default=0, metavar="INT",
                              help="Keep original quality scores in output (0/1) [0]")
    output_group.add_argument("--phred_alpha", type=str, default='phred42', metavar="STR",
                              help="Phred quality (q-score) ascii character alphabet used by input (phred42, phred63, phred94) [phred42]")

    # Performance Group
    perf_group = parser.add_argument_group("PERFORMANCE & PARALLELIZATION")
    perf_group.add_argument("--workers", type=int, default=1, metavar="INT",
                            help="Number of parallel workers (use 4+ for large files >5GB) [1]")
    perf_group.add_argument("--chunk_mb", type=int, default=32, metavar="INT",
                            help="Chunk size in MB for parallel processing [32]")
    perf_group.add_argument("--profile", type=int, default=0, metavar="INT",
                            help="Enable profiling (0/1) [0]")

    args = parser.parse_args()
    
    
    if args.qual_scale == 'custom' and args.custom_formula is None:
        print("ERROR: --qual_scale custom requires --custom_formula argument")
        print("\nExample usage:")
        print("  --qual_scale custom --custom_formula '1 + 62 * (x - 40) / 53'")
        print("  --qual_scale custom --custom_formula 'ln(x - 39) / ln(54) * 62 + 1'")
        exit(1)
    
    if args.multi_flow == 1 and args.compress_hdr == 0:
        print("WARNING: --multi_flow requires --compress_hdr 1 to function")
        print("Enabling header compression automatically...")
        args.compress_hdr = 1
    
    if args.rm_repeat_hdr == 1 and args.compress_hdr == 0:
        print("WARNING: --rm_repeat_hdr requires --compress_hdr 1 to function")
        print("Enabling header compression automatically...")
        args.compress_hdr = 1


    # Mode shortcuts
    if args.mode == 0:
        args.compress_hdr = 1
        args.bin_write = 0
        args.keep_bases = 1
        args.keep_qual = 1
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
            args.custom_formula,
            phred_alphabet_max
        )

    # Create base map
    numpy_base_map = np.zeros(128, dtype=np.uint8)
    numpy_base_map[ord("N")] = args.gray_N
    numpy_base_map[ord("A")] = args.gray_A
    numpy_base_map[ord("T")] = args.gray_T
    numpy_base_map[ord("C")] = args.gray_C
    numpy_base_map[ord("G")] = args.gray_G

    start_time = time.perf_counter()

    # Initialize profiler if requested
    profiler = None
    if args.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        print("Profiling enabled...")

    # Create phred map if needed
    phred_map = create_phred_quality_map(args.phred_off, phred_alphabet_max) if args.extract_qual else None

    print(f"Converting sequences from {args.input_path}...")
    
    export_scalars_to_txt(
        args.input_path, numpy_base_map,
        args.output_path, phred_map, args.min_qual,
        args.qual_scale, args.bin_write,
        compress_headers=(args.compress_hdr == 1),
        sequencer_type=args.seq_type,
        sra_accession=args.sra_acc,
        keep_bases=args.keep_bases, keep_quality=args.keep_qual,
        custom_formula=args.custom_formula,
        multiple_flowcells=(args.multi_flow == 1),
        remove_repeating_header=(args.rm_repeat_hdr == 1),
        phred_alphabet_max=phred_alphabet_max,
        paired_end=(args.paired == 1),
        paired_end_mode=args.paired_mode,
        num_workers=args.workers,
        chunk_size_mb=args.chunk_mb
    )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Task completed in {elapsed_time:.4f} seconds")
    print(f"Saved scalars to {args.output_path}")
    
    if profiler is not None:
        profiler.disable()
        print("\n" + "="*80)
        print("Profiling Results:")
        print("="*80)
        stats = pstats.Stats(profiler)
        stats.sort_stats('cumulative')
        stats.print_stats(30)
        stats.sort_stats('tottime')
        print("\n" + "="*80)
        print("By total time:")
        print("="*80)
        stats.print_stats(30)


if __name__ == "__main__":
    main()
