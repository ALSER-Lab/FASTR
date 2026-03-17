#!/bin/bash
set -e

# Change to performance directory
cd "$(dirname "$0")"

echo "=========================================="
echo "FASTR Format Test - Fixed Version"
echo "=========================================="
echo ""

# First, decompress test files if needed
echo "Preparing test data..."
if [ ! -f test_data/sample2_1_trimmed_1.fastq ]; then
    echo "Decompressing R1..."
    zcat test_data/sample2_1_trimmed_1.fastq.gz > test_data/sample2_1_trimmed_1.fastq
fi
if [ ! -f test_data/sample2_1_trimmed_2.fastq ]; then
    echo "Decompressing R2..."
    zcat test_data/sample2_1_trimmed_2.fastq.gz > test_data/sample2_1_trimmed_2.fastq
fi

# Create FASTR files
echo "Step 1: Converting FASTQ to FASTR (Mode 1)..."
time python ../../src/toFASTR.py \
    test_data/sample2_1_trimmed_1.fastq \
    test_data/sample2_1_trimmed_1_fixed.fastr \
    --mode 1 \
    --threads 1 2>&1 | grep -E "(Total sequences|completed)"

time python ../../src/toFASTR.py \
    test_data/sample2_1_trimmed_2.fastq \
    test_data/sample2_1_trimmed_2_fixed.fastr \
    --mode 1 \
    --threads 1 2>&1 | grep -E "(Total sequences|completed)"

echo ""
echo "Step 2: Streaming FASTR → FASTQ → BWA-MEM2..."
time (bwa-mem2 mem -t 4 test_data/chr22.fasta \
    <(python ../../src/toFASTQ.py test_data/sample2_1_trimmed_1_fixed.fastr /dev/stdout --threads 1 2>/dev/null) \
    <(python ../../src/toFASTQ.py test_data/sample2_1_trimmed_2_fixed.fastr /dev/stdout --threads 1 2>/dev/null) \
    | samtools view -b -o results/test_fastr_fixed_streaming.bam) 2>&1 | grep -A 20 "Overall time"

echo ""
echo "Step 3: Checking alignment results..."
samtools flagstat results/test_fastr_fixed_streaming.bam

echo ""
echo "Step 4: Comparing with baseline..."
echo "Baseline (internal reader):  545,788 alignments"
echo "Pigz streaming:             545,788 alignments"
echo "Fixed FASTR streaming:      $(samtools view -c results/test_fastr_fixed_streaming.bam) alignments"

echo ""
echo "=========================================="
echo "FASTR Test Complete!"
echo "=========================================="
