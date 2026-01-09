#!/bin/bash
# Usage: ./launch_benchmark.sh /path/to/input.fastq /path/to/output_dir

INPUT_FASTQ="$1"
THREADS=64
OUTDIR="$2"

if [[ -z "$INPUT_FASTQ" ]]; then
    echo "Usage: ./launch_benchmark.sh /path/to/input.fastq /path/to/output_dir"
    exit 1
fi

mkdir -p "$OUTDIR"
BASENAME=$(basename "$INPUT_FASTQ" .fastq)
BASENAME=$(basename "$BASENAME" .fq)

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
echo "All jobs submitted. Check status with: squeue -u \$USER"
