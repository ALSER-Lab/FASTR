#!/bin/bash
# Setup script for test data
# Decompresses test files and builds BWA index if needed

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/test_data"

echo "=========================================="
echo "Setting up FASTR test data"
echo "=========================================="
echo ""

# Check if compressed files exist
if [ ! -f "$DATA_DIR/sample2_1_trimmed_1.fastq.gz" ]; then
    echo "❌ Error: test_data/sample2_1_trimmed_1.fastq.gz not found"
    exit 1
fi

if [ ! -f "$DATA_DIR/chr22.fasta.gz" ]; then
    echo "❌ Error: test_data/chr22.fasta.gz not found"
    exit 1
fi

# Decompress FASTQ files if not already done
echo "Step 1: Decompressing FASTQ files..."
if [ ! -f "$DATA_DIR/sample2_1_trimmed_1.fastq" ]; then
    echo "  Decompressing R1..."
    zcat "$DATA_DIR/sample2_1_trimmed_1.fastq.gz" > "$DATA_DIR/sample2_1_trimmed_1.fastq"
else
    echo "  R1 already decompressed"
fi

if [ ! -f "$DATA_DIR/sample2_1_trimmed_2.fastq" ]; then
    echo "  Decompressing R2..."
    zcat "$DATA_DIR/sample2_1_trimmed_2.fastq.gz" > "$DATA_DIR/sample2_1_trimmed_2.fastq"
else
    echo "  R2 already decompressed"
fi

# Decompress reference genome if not already done
echo ""
echo "Step 2: Decompressing reference genome..."
if [ ! -f "$DATA_DIR/chr22.fasta" ]; then
    echo "  Decompressing chr22.fasta..."
    zcat "$DATA_DIR/chr22.fasta.gz" > "$DATA_DIR/chr22.fasta"
else
    echo "  chr22.fasta already decompressed"
fi

# Build BWA index if not already done
echo ""
echo "Step 3: Building BWA-MEM2 index (if needed)..."
if [ ! -f "$DATA_DIR/chr22.fasta.bwt.2bit.64" ]; then
    echo "  Building index... (this may take 1-2 minutes)"
    cd "$DATA_DIR"
    bwa-mem2 index chr22.fasta
    echo "  Index built successfully"
else
    echo "  Index already exists"
fi

echo ""
echo "=========================================="
echo "Test data setup complete!"
echo "=========================================="
echo ""
echo "Available files:"
ls -lh "$DATA_DIR"/*.fastq "$DATA_DIR"/chr22.fasta 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""
