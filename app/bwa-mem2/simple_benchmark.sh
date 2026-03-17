#!/bin/bash -euo pipefail

# Simple benchmark script to compare bwa-mem2 performance
# Run this with: pixi run benchmark-simple

echo "==================================================================="
echo "BWA-MEM2 Performance Benchmark"
echo "==================================================================="

# Change to performance directory
cd "$(dirname "$0")"

# Find the BWA-MEM2 index prefix
INDEX=$(find -L test_data/ -name "*.amb" | sed 's/\.amb$//')
echo "Index: $INDEX"

# Input files
READ1="test_data/sample2_1_trimmed_1.fastq.gz"
READ2="test_data/sample2_1_trimmed_2.fastq.gz"

# Read group
RG="@RG\tID:sample2_1\tSM:sample2\tPL:ILLUMINA\tLB:sample2_lib\tPU:1"

# Number of threads
THREADS=4

echo ""
echo "Input files:"
ls -lh $READ1 $READ2

echo ""
echo "==================================================================="
echo "Test 1: BWA-MEM2 with internal reader (current method)"
echo "==================================================================="
echo "Command: bwa-mem2 mem -t $THREADS $INDEX $READ1 $READ2 | samtools view -Sb > test1_internal.bam"
echo ""

time (bwa-mem2 mem \
    -t $THREADS \
    -M \
    -R "$RG" \
    $INDEX \
    $READ1 $READ2 \
    | samtools view -Sb - > test1_internal.bam)

echo ""
echo "Output: $(ls -lh test1_internal.bam | awk '{print $5, $9}')"

echo ""
echo "==================================================================="
echo "Test 2: BWA-MEM2 with pigz streaming (process substitution)"
echo "==================================================================="
echo "Command: bwa-mem2 mem -t $THREADS $INDEX <(pigz -cd -p 2 $READ1) <(pigz -cd -p 2 $READ2) | samtools view -Sb > test2_pigz.bam"
echo ""

time (bwa-mem2 mem \
    -t $THREADS \
    -M \
    -R "$RG" \
    $INDEX \
    <(pigz -cd -p 2 $READ1) \
    <(pigz -cd -p 2 $READ2) \
    | samtools view -Sb - > test2_pigz.bam)

echo ""
echo "Output: $(ls -lh test2_pigz.bam | awk '{print $5, $9}')"

echo ""
echo "==================================================================="
echo "Verification"
echo "==================================================================="

# Quick check that both produce same number of alignments
COUNT1=$(samtools view -c test1_internal.bam)
COUNT2=$(samtools view -c test2_pigz.bam)

echo "Test 1 alignments: $COUNT1"
echo "Test 2 alignments: $COUNT2"

if [ "$COUNT1" -eq "$COUNT2" ]; then
    echo "✓ Both methods produced the same number of alignments"
else
    echo "✗ Warning: Different number of alignments!"
fi

echo ""
echo "==================================================================="
echo "Summary"
echo "==================================================================="
echo "Review the timing information above to determine which is faster."
echo ""
echo "Cleanup test files with: rm -f test1_internal.bam test2_pigz.bam"
