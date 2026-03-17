#!/bin/bash
set -e

# Change to performance directory
cd "$(dirname "$0")"

echo "=== Verifying Original FASTQ Files Are Properly Paired ==="
echo ""

echo "Step 1: Count total reads in each file"
R1_COUNT=$(zcat test_data/sample2_1_trimmed_1.fastq.gz | wc -l | awk '{print $1/4}')
R2_COUNT=$(zcat test_data/sample2_1_trimmed_2.fastq.gz | wc -l | awk '{print $1/4}')

echo "R1 reads: $R1_COUNT"
echo "R2 reads: $R2_COUNT"

if [ "$R1_COUNT" != "$R2_COUNT" ]; then
    echo "❌ ERROR: Read counts don't match!"
    exit 1
else
    echo "✅ Read counts match"
fi

echo ""
echo "Step 2: Check if read IDs match between R1 and R2"
echo "Extracting first 100 read IDs from R1..."
zcat test_data/sample2_1_trimmed_1.fastq.gz | awk 'NR % 4 == 1' | head -100 > /tmp/r1_ids.txt

echo "Extracting first 100 read IDs from R2..."
zcat test_data/sample2_1_trimmed_2.fastq.gz | awk 'NR % 4 == 1' | head -100 > /tmp/r2_ids.txt

echo ""
echo "Comparing read IDs (should be identical)..."
if diff /tmp/r1_ids.txt /tmp/r2_ids.txt > /dev/null 2>&1; then
    echo "✅ First 100 read IDs match perfectly!"
else
    echo "⚠️ Read IDs differ. Let's check the differences:"
    diff /tmp/r1_ids.txt /tmp/r2_ids.txt | head -20
fi

echo ""
echo "Step 3: Sample some read IDs from different positions"
echo "First 3 pairs:"
echo "R1:" && zcat test_data/sample2_1_trimmed_1.fastq.gz | awk 'NR % 4 == 1' | head -3
echo "R2:" && zcat test_data/sample2_1_trimmed_2.fastq.gz | awk 'NR % 4 == 1' | head -3

echo ""
echo "Middle 3 pairs (around read 136000):"
echo "R1:" && zcat test_data/sample2_1_trimmed_1.fastq.gz | awk 'NR % 4 == 1' | sed -n '136000,136002p'
echo "R2:" && zcat test_data/sample2_1_trimmed_2.fastq.gz | awk 'NR % 4 == 1' | sed -n '136000,136002p'

echo ""
echo "Last 3 pairs:"
echo "R1:" && zcat test_data/sample2_1_trimmed_1.fastq.gz | awk 'NR % 4 == 1' | tail -3
echo "R2:" && zcat test_data/sample2_1_trimmed_2.fastq.gz | awk 'NR % 4 == 1' | tail -3

echo ""
echo "Step 4: Verify with BWA-MEM2 (should work fine with original files)"
echo "Running quick alignment test with original files..."
bwa-mem2 mem -t 4 test_data/chr22.fasta \
    <(zcat test_data/sample2_1_trimmed_1.fastq.gz | head -4000) \
    <(zcat test_data/sample2_1_trimmed_2.fastq.gz | head -4000) \
    2>&1 | grep -E "(read [0-9]+ sequences|paired reads have different names|analyzing insert size)" | head -10

echo ""
echo "=== Verification Complete ==="
