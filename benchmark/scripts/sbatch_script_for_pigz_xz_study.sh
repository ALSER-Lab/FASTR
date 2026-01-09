#!/bin/bash
#===============================================================================
# FASTQ Compression Benchmark Script
#===============================================================================
#
# DESCRIPTION:
#   Benchmarks compression and decompression performance of pigz (parallel gzip)
#   and xz on FASTQ files using a SLURM HPC cluster. Compares speed, compression
#   ratio, and resource usage across different compression levels.
#
# COMPRESSION TOOLS TESTED:
#   - pigz -1  (fast)    : Fastest gzip, lowest compression
#   - pigz -9  (best)    : Best gzip compression, slower
#   - pigz -11 (zopfli)  : Maximum compression using Zopfli algorithm (very slow)
#   - xz -1    (fast)    : Fast LZMA compression
#   - xz -9    (best)    : Best LZMA compression (slow but smallest files)
#
# WORKFLOW:
#   1. Submits 5 compression jobs in parallel
#   2. Each decompression job waits for its corresponding compression to finish
#   3. All jobs use /usr/bin/time -v to capture detailed resource metrics
#
# OUTPUT FILES:
#   Compressed files:
#     - <basename>.pigz_fast.gz      (pigz -1)
#     - <basename>.pigz_best.gz      (pigz -9)
#     - <basename>.pigz_zopfli.gz    (pigz -11)
#     - <basename>.xz_fast.xz        (xz -1)
#     - <basename>.xz_best.xz        (xz -9)
#
#   Decompressed files:
#     - <basename>.<method>.decompressed.fastq
#
#   Log files (contain timing and memory stats):
#     - <method>_compress_<jobid>.out/.err
#     - <method>_decompress_<jobid>.out/.err
#
# SLURM SETTINGS:
#   - Partition: qDEV
#   - CPUs: 64 threads per job
#   - Time limit: 11 hours per job
#
# USAGE:
#   ./run_pigz_xz.sh <input.fastq> <output_directory>
#
# EXAMPLE:
#   ./run_pigz_xz.sh /data/sample.fastq /results/compression_benchmark
#
# MONITORING:
#   squeue -u $USER                    # Check job status
#   tail -f <outdir>/*_compress_*.out  # Watch compression progress
#
# PARSING RESULTS:
#   After completion, compare file sizes:
#     ls -lh <outdir>/*.gz <outdir>/*.xz
#
#   Extract timing from logs:
#     grep "Elapsed" <outdir>/*.out
#
#   Extract peak memory:
#     grep "Maximum resident" <outdir>/*.out
#
# DEPENDENCIES:
#   - pigz (parallel gzip)
#   - xz (with threading support)
#   - SLURM workload manager
#===============================================================================

INPUT_FASTQ="$1"
THREADS=64
OUTDIR="$2"

if [[ -z "$INPUT_FASTQ" || -z "$OUTDIR" ]]; then
    echo "Usage: ./run_pigz_xz.sh /path/to/input.fastq /path/to/output_dir"
    exit 1
fi

if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file not found: $INPUT_FASTQ"
    exit 1
fi

mkdir -p "$OUTDIR"
BASENAME=$(basename "$INPUT_FASTQ" .fastq)
BASENAME=$(basename "$BASENAME" .fq)

echo "=============================================="
echo "FASTQ Compression Benchmark"
echo "=============================================="
echo "Input:   $INPUT_FASTQ"
echo "Output:  $OUTDIR"
echo "Threads: $THREADS"
echo "=============================================="
echo ""

# Submit compression jobs
JOB1=$(sbatch --parsable --job-name=xz_fast_c --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/xz_fast_compress_%j.out --error=${OUTDIR}/xz_fast_compress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'xz -1 -k -T $THREADS -c $INPUT_FASTQ > ${OUTDIR}/${BASENAME}.xz_fast.xz' 2>&1")
echo "Submitted xz_fast compress: $JOB1"

JOB2=$(sbatch --parsable --job-name=xz_best_c --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/xz_best_compress_%j.out --error=${OUTDIR}/xz_best_compress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'xz -9 -k -T $THREADS -c $INPUT_FASTQ > ${OUTDIR}/${BASENAME}.xz_best.xz' 2>&1")
echo "Submitted xz_best compress: $JOB2"

JOB3=$(sbatch --parsable --job-name=pigz_fast_c --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_fast_compress_%j.out --error=${OUTDIR}/pigz_fast_compress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -1 -k -p $THREADS -c $INPUT_FASTQ > ${OUTDIR}/${BASENAME}.pigz_fast.gz' 2>&1")
echo "Submitted pigz_fast compress: $JOB3"

JOB4=$(sbatch --parsable --job-name=pigz_best_c --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_best_compress_%j.out --error=${OUTDIR}/pigz_best_compress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -9 -k -p $THREADS -c $INPUT_FASTQ > ${OUTDIR}/${BASENAME}.pigz_best.gz' 2>&1")
echo "Submitted pigz_best compress: $JOB4"

JOB5=$(sbatch --parsable --job-name=pigz_zopfli_c --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_zopfli_compress_%j.out --error=${OUTDIR}/pigz_zopfli_compress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -11 -k -p $THREADS -c $INPUT_FASTQ > ${OUTDIR}/${BASENAME}.pigz_zopfli.gz' 2>&1")
echo "Submitted pigz_zopfli compress: $JOB5"

# Submit decompression jobs with dependency on compression
sbatch --dependency=afterok:$JOB1 --job-name=xz_fast_d --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/xz_fast_decompress_%j.out --error=${OUTDIR}/xz_fast_decompress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'xz -d -k -T $THREADS -c ${OUTDIR}/${BASENAME}.xz_fast.xz > ${OUTDIR}/${BASENAME}.xz_fast.decompressed.fastq' 2>&1"
echo "Submitted xz_fast decompress (depends on $JOB1)"

sbatch --dependency=afterok:$JOB2 --job-name=xz_best_d --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/xz_best_decompress_%j.out --error=${OUTDIR}/xz_best_decompress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'xz -d -k -T $THREADS -c ${OUTDIR}/${BASENAME}.xz_best.xz > ${OUTDIR}/${BASENAME}.xz_best.decompressed.fastq' 2>&1"
echo "Submitted xz_best decompress (depends on $JOB2)"

sbatch --dependency=afterok:$JOB3 --job-name=pigz_fast_d --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_fast_decompress_%j.out --error=${OUTDIR}/pigz_fast_decompress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -d -k -p $THREADS -c ${OUTDIR}/${BASENAME}.pigz_fast.gz > ${OUTDIR}/${BASENAME}.pigz_fast.decompressed.fastq' 2>&1"
echo "Submitted pigz_fast decompress (depends on $JOB3)"

sbatch --dependency=afterok:$JOB4 --job-name=pigz_best_d --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_best_decompress_%j.out --error=${OUTDIR}/pigz_best_decompress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -d -k -p $THREADS -c ${OUTDIR}/${BASENAME}.pigz_best.gz > ${OUTDIR}/${BASENAME}.pigz_best.decompressed.fastq' 2>&1"
echo "Submitted pigz_best decompress (depends on $JOB4)"

sbatch --dependency=afterok:$JOB5 --job-name=pigz_zopfli_d --partition=qDEV --ntasks=1 --cpus-per-task=64 --time=11:00:00 \
    --output=${OUTDIR}/pigz_zopfli_decompress_%j.out --error=${OUTDIR}/pigz_zopfli_decompress_%j.err \
    --wrap="/usr/bin/time -v sh -c 'pigz -d -k -p $THREADS -c ${OUTDIR}/${BASENAME}.pigz_zopfli.gz > ${OUTDIR}/${BASENAME}.pigz_zopfli.decompressed.fastq' 2>&1"
echo "Submitted pigz_zopfli decompress (depends on $JOB5)"

echo ""
echo "=============================================="
echo "All jobs submitted!"
echo "=============================================="
echo "Monitor with: squeue -u \$USER"
echo "Results in:   $OUTDIR"
echo "=============================================="
