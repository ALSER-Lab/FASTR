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
                          quality_scaling='none', binary=True, log_a=None, compress_headers=False,
                          sequencer_type='none', sra_accession=None, keep_bases=False, keep_quality=False,
                          binary_bases=False, binary_quality=False, custom_formula=None, multiple_flowcells=False,
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
        if compress_headers and sequencer_type != 'none':
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
                binary_bases=binary_bases,
                keep_quality=keep_quality,
                binary_quality=binary_quality,
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
        if compress_headers and sequencer_type != 'none':
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
    argument_parser = argparse.ArgumentParser(description="Convert and compress FASTQ/FASTA files to scalar format")
    argument_parser.add_argument("input_path", type=str, help="Path of .fasta or .fastq file")
    argument_parser.add_argument("output_path", type=str, help="Output file path")

    argument_parser.add_argument("--quality_scaling", type=str,
                                 choices=['log', 'log_reverse', 'log_custom', 'custom'],
                                 default='none',
                                 help="Quality scaling method (default: none)")

    # Quick mode implementation
    argument_parser.add_argument("--mode", type=int,
                            help="What mode to use when writing out converted values. "
                            "1: Header compression only "
                            "2: Base conversion into numbers only "
                            "3: Header and base conversion, written out in two lines "
                            "4: Repeating header removal entirely, base conversion kept, written out in one line")
    
    # Paired-end args
    argument_parser.add_argument("--paired_end", type=int, default=0,
                                help="Signifies if the input file contains paired-end reads (0/1, default 0)")
    argument_parser.add_argument("--paired_end_mode", type=str,
                                 choices=['same_file', 'separate_files'],
                                 default='same_file',
                                 help="Output mode for paired-end reads (default: same_file)")
    
    # SRA accession
    argument_parser.add_argument("--sra_accession", type=str, default=None,
                                help="SRA accession number (e.g., SRR12345678)")
    
    # Header compression arguments
    argument_parser.add_argument("--compress_headers", type=int, default=0,
                                help="Compress FASTQ headers on-the-fly (0/1, default 0)")
    argument_parser.add_argument("--sequencer_type", type=str,
                                choices=['none', 'illumina', 'pacbio_ccs', 'pacbio_hifi', 'pacbio_subread', 'pacbio_clr', 'ont', 'srr'],
                                default='none',
                                help="Sequencer type for header compression (default: none)")
    argument_parser.add_argument("--multiple_flowcells", type=int, default=0,
                                help="Enable multiple flowcell detection and tracking (0/1, default 0)")
    
    # Quality arguments
    argument_parser.add_argument("--extract_quality", type=int, default=1,
                                help="For FASTQ: extract quality scores (0/1, default 1)")
    argument_parser.add_argument("--phred_offset", type=int, default=33,
                                help="Phred quality offset (default 33)")
    argument_parser.add_argument("--min_quality", type=int, default=0,
                                help="Minimum quality score threshold (default 0)")
    argument_parser.add_argument("--custom_formula", type=str, default=None,
                                help="Custom formula for quality scaling (use 'x' for quality score). "
                                     "Example: '1 + 62 * (x - 40) / 53' or 'ln(x) * 10'")
    argument_parser.add_argument("--log_a", type=int, default=0,
                                help="Tunable a value for application in custom log function (default 0)")
    
    # Output format
    argument_parser.add_argument("--binary_write", type=int, default=1,
                                help="Enable binary writing of sequence integers (0/1, default 1)")
    argument_parser.add_argument("--keep_bases", type=int, default=0,
                                help="Return textual bases without scaling or one-hot encoding (0/1, default 0)")
    argument_parser.add_argument("--keep_quality", type=int, default=0,
                            help="Keep original quality scores in output (0/1, default 0)")
    argument_parser.add_argument("--binary_bases", type=int, default=0,
                            help="Use binary encoding for STRING BASES (0/1, default 0)")
    argument_parser.add_argument("--binary_quality", type=int, default=0,
                            help="Write quality scores as binary numeric values (0/1, default 0)")
    
    argument_parser.add_argument("--remove_repeating_header", type=int, default=0,
                            help="Remove repeating metadata from headers, store only at top (0/1, default 0)")
    argument_parser.add_argument("--phred_alphabet", type=str, default='phred42',
                            help="Phred quality (q-score) ascii character alphabet used by the input fastq")

    # Base mapping
    argument_parser.add_argument("--gray_N", type=int, default=1, help="Grayscale value for N (default 1)")
    argument_parser.add_argument("--gray_A", type=int, default=63, help="Grayscale value for A (default 63)")
    argument_parser.add_argument("--gray_T", type=int, default=127, help="Grayscale value for T (default 127)")
    argument_parser.add_argument("--gray_C", type=int, default=191, help="Grayscale value for C (default 191)")
    argument_parser.add_argument("--gray_G", type=int, default=255, help="Grayscale value for G (default 255)")
    
    # Profiling
    argument_parser.add_argument("--profile", type=int, default=0,
                                help="Enable profiling (0/1, default 0)")
    
    # Multiprocessing
    argument_parser.add_argument("--num_workers", type=int, default=1,
                                help="Number of parallel workers (default: 1, use 4+ for large files >5GB)")
    
    argument_parser.add_argument("--chunk_size_mb", type=int, default=32,
                                help="Chunk size in MB for parallel processing (default: 32)")

    user_arguments = argument_parser.parse_args()
    
    # Validation
    if user_arguments.keep_bases == 1 and user_arguments.binary_bases == 1:
        print("ERROR: Cannot use both --keep_bases and --binary_bases. Choose one:")
        print("  --keep_bases 1: Keep ASCII letters (A=65, T=84, etc.)")
        print("  --binary_bases 1: Convert to compact encoding (A=0, T=1, etc.)")
        exit(1)
    
    if user_arguments.quality_scaling == 'custom' and user_arguments.custom_formula is None:
        print("ERROR: --quality_scaling custom requires --custom_formula argument")
        print("\nExample usage:")
        print("  --quality_scaling custom --custom_formula '1 + 62 * (x - 40) / 53'")
        print("  --quality_scaling custom --custom_formula 'ln(x - 39) / ln(54) * 62 + 1'")
        exit(1)
    
    if user_arguments.multiple_flowcells == 1 and user_arguments.compress_headers == 0:
        print("WARNING: --multiple_flowcells requires --compress_headers 1 to function")
        print("Enabling header compression automatically...")
        user_arguments.compress_headers = 1
    
    if user_arguments.remove_repeating_header == 1 and user_arguments.compress_headers == 0:
        print("WARNING: --remove_repeating_header requires --compress_headers 1 to function")
        print("Enabling header compression automatically...")
        user_arguments.compress_headers = 1

    # Mode shortcuts
    if user_arguments.mode == 1:
        user_arguments.compress_headers = 1
        user_arguments.binary_write = 0
    elif user_arguments.mode == 2:
        user_arguments.compress_headers = 0
        user_arguments.binary_write = 1
    elif user_arguments.mode == 3:
        user_arguments.compress_headers = 1
        user_arguments.binary_write = 1
    elif user_arguments.mode == 4:
        user_arguments.compress_headers = 1
        user_arguments.remove_repeating_header = 1
        user_arguments.binary_write = 1

    # Phred alphabet configuration
    if user_arguments.phred_alphabet == "phred42":
        phred_alphabet_max = 41
    elif user_arguments.phred_alphabet == "phred63":
        phred_alphabet_max = 62
    elif user_arguments.phred_alphabet == "phred94":
        phred_alphabet_max = 93

    # Validate custom formula if provided
    if user_arguments.custom_formula:
        user_arguments.custom_formula = validate_and_adjust_formula(
            user_arguments.custom_formula,
            phred_alphabet_max
        )

    # Create base map
    numpy_base_map = np.zeros(128, dtype=np.uint8)
    numpy_base_map[ord("N")] = user_arguments.gray_N
    numpy_base_map[ord("A")] = user_arguments.gray_A
    numpy_base_map[ord("T")] = user_arguments.gray_T
    numpy_base_map[ord("C")] = user_arguments.gray_C
    numpy_base_map[ord("G")] = user_arguments.gray_G

    start_time = time.perf_counter()

    # Initialize profiler if requested
    profiler = None
    if user_arguments.profile == 1:
        profiler = cProfile.Profile()
        profiler.enable()
        print("Profiling enabled...")

    # Create phred map if needed
    phred_map = create_phred_quality_map(user_arguments.phred_offset, phred_alphabet_max) if user_arguments.extract_quality else None

    print(f"Converting sequences from {user_arguments.input_path}...")
    
    export_scalars_to_txt(
        user_arguments.input_path, numpy_base_map,
        user_arguments.output_path, phred_map, user_arguments.min_quality,
        user_arguments.quality_scaling, user_arguments.binary_write,
        user_arguments.log_a, compress_headers=(user_arguments.compress_headers == 1),
        sequencer_type=user_arguments.sequencer_type,
        sra_accession=user_arguments.sra_accession,
        keep_bases=user_arguments.keep_bases, keep_quality=user_arguments.keep_quality,
        binary_bases=user_arguments.binary_bases, binary_quality=user_arguments.binary_quality,
        custom_formula=user_arguments.custom_formula,
        multiple_flowcells=(user_arguments.multiple_flowcells == 1),
        remove_repeating_header=(user_arguments.remove_repeating_header == 1),
        phred_alphabet_max=phred_alphabet_max,
        paired_end=(user_arguments.paired_end == 1),
        paired_end_mode=user_arguments.paired_end_mode,
        num_workers=user_arguments.num_workers,
        chunk_size_mb=user_arguments.chunk_size_mb
    )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Task completed in {elapsed_time:.4f} seconds")
    print(f"Saved scalars to {user_arguments.output_path}")
    
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