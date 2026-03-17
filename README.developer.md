# FASTR Performance Testing - Developer Guide

This directory contains comprehensive performance benchmarking and testing setup for the FASTR (Fast Aligned Sequencing Transformation and Reconstruction) format, including:

- Bug fix verification for critical data loss issue in FASTQ to FASTR conversion
- Performance comparison of different BWA-MEM2 alignment methods
- Complete testing setup with reproducible environment using Pixi

## Quick Start

### 1. Install Pixi (one-time setup)

```bash
curl -fsSL https://pixi.sh/install.sh | bash
# Add pixi to PATH (follow installer instructions)
```

### 2. Clone and Setup FASTR with Performance Testing

```bash
git clone https://github.com/ALSER-Lab/FASTR.git
cd FASTR
```

### 3. Initialize Pixi Environment

```bash
# Create and activate the pixi environment
pixi install

# Install FASTR in development mode
pixi run install-dev

# Setup test data (decompress, build BWA index)
pixi run setup-test-data
```

### 4. Verify Everything Works

```bash
# Verify test data is set up
ls -lh app/bwa-mem2/test_data/

# Check that all dependencies are available
pixi run verify-pairing
```

## Test Data

The `app/bwa-mem2/test_data/` directory contains:

| File | Size | Type | Notes |
|------|------|------|-------|
| `sample2_1_trimmed_1.fastq.gz` | 29 MB | Compressed | R1 reads (272,479 sequences) |
| `sample2_1_trimmed_2.fastq.gz` | 32 MB | Compressed | R2 reads (272,479 sequences) |
| `chr22.fasta.gz` | 11 MB | Compressed | Reference genome (chr22) |
| **Total** | **72 MB** | | Fits in GitHub Actions cache |

**First-time setup:** The `setup-test-data` task will:
1. Decompress FASTQ files to working directory
2. Decompress reference genome
3. Build BWA-MEM2 index (takes ~1-2 minutes, only once)

**Subsequent runs:** Index and decompressed files are reused for faster execution.

## Available Tasks

All tasks can be run with `pixi run <task-name>`:

### `pixi run setup-test-data`
Prepares test data for benchmarking and testing.

**What it does:**
- Decompresses FASTQ files from .gz archives
- Decompresses reference genome
- Builds BWA-MEM2 index for the reference genome
- ~2 minutes execution time (first run only)

**Note:** Only needs to be run once. Subsequent runs detect existing files and skip decompression/indexing.

**Output:**
```
Step 1: Decompressing FASTQ files...
Step 2: Decompressing reference genome...
Step 3: Building BWA-MEM2 index (if needed)...
```

### `pixi run verify-pairing`
Verifies that the input FASTQ files are properly paired and synchronized.

**What it does:**
- Counts reads in both R1 and R2 files
- Verifies read IDs match at beginning, middle, and end
- Quick test with BWA-MEM2 to confirm proper pairing
- ~1 minute execution time

**Output:**
```
R1 reads: 272479
R2 reads: 272479
✅ Read counts match
✅ First 100 read IDs match perfectly!
```

### `pixi run benchmark-simple`
Compares basic BWA-MEM2 performance between internal reader and pigz streaming.

**What it does:**
- Method 1: BWA-MEM2 with internal reader (baseline)
- Method 2: BWA-MEM2 with pigz multi-threaded decompression
- Compares timing and output file sizes
- ~3 minutes execution time

**Output:**
```
Test 1 alignments: 545788
Test 2 alignments: 545788
✓ Both methods produced the same number of alignments
```

### `pixi run benchmark-fastr`
Tests the fixed FASTR format conversion and streaming with BWA-MEM2.

**What it does:**
- Converts FASTQ to FASTR (Mode 1)
- Streams FASTR back to FASTQ
- Runs alignment via piped streaming
- Verifies all reads are preserved (272,479 reads)
- ~5 minutes execution time

**Output:**
```
Total sequences written: 272,479
545788 + 0 in total (QC-passed reads)
Fixed FASTR streaming: 545788 alignments
```

### `pixi run benchmark-all`
Runs both benchmark-simple and benchmark-fastr sequentially.

### `pixi run install-dev`
Installs FASTR package in development mode for local testing.

## Directory Structure

```
FASTR/
├── pixi.toml                 # Pixi environment configuration
├── README.developer.md       # This file
├── app/
│   ├── bwa-mem2/
│   │   ├── test_data/           # Test datasets (FASTQ, reference genome, indexes)
│   │   │   ├── sample2_1_trimmed_1.fastq.gz      # Test read file 1 (29MB)
│   │   │   ├── sample2_1_trimmed_2.fastq.gz      # Test read file 2 (32MB)
│   │   │   ├── chr22.fasta                       # Reference genome (49MB)
│   │   │   └── chr22.fasta.*                     # BWA-MEM2 index files
│   │   ├── results/              # Output files from benchmark runs
│   │   │   ├── test1_internal.bam
│   │   │   ├── test2_pigz.bam
│   │   │   └── test_fastr_fixed_streaming.bam
│   │   ├── docs/                # Documentation
│   │   │   ├── FASTR_BUG_FIX_SUMMARY.md           # Bug details and fix
│   │   │   ├── FINAL_RESULTS_WITH_FIXED_FASTR.md  # Complete results
│   │   │   └── BENCHMARK_RESULTS.md               # Initial benchmarks
│   │   ├── simple_benchmark.sh   # 2-method performance comparison
│   │   ├── test_fastr_fixed.sh   # FASTR streaming test
│   │   ├── verify_pairing.sh     # FASTQ validation
│   │   └── setup_test_data.sh    # Test data setup
│   └── ...other apps...
├── src/
│   ├── toFASTR.py           # FASTQ → FASTR converter
│   ├── toFASTQ.py           # FASTR → FASTQ reconstructor
│   └── toFASTR_chunk_processor.py  # [FIXED] Chunk boundary detection
└── ...other files...
```

## The FASTR Bug Fix

### What Was the Bug?

The original `toFASTR_chunk_processor.py` used an **unreliable method** to find FASTQ record boundaries:

```python
# BUGGY CODE (line 135):
last_at = buffer.rfind(b"\n@")
```

**Problem:** Quality scores in FASTQ can contain "@" characters (ASCII 64). When a quality line has `\n@`, the code incorrectly splits the chunk, dropping reads.

**Impact:**
- Input: 272,479 reads
- Output: 272,472 reads (7 reads lost!)
- Files become desynchronized
- BWA-MEM2 produces 0 alignments

### How Was It Fixed?

Implemented proper FASTQ record structure validation:

```python
def find_last_fastq_record_boundary(buffer: bytes) -> int:
    """
    Find position of last complete FASTQ record by validating:
    1. Header line starts with @
    2. Separator line starts with +
    3. Quality length matches sequence length
    """
    # ... implementation ...
```

**Result After Fix:**
- Input: 272,479 reads
- Output: 272,479 reads ✅
- All reads preserved with correct pairing
- BWA-MEM2 produces 545,788 alignments ✅

## Test Data

The benchmark uses real human reads from **chr22** with proper pairing:

| File | Size | Reads | Type |
|------|------|-------|------|
| sample2_1_trimmed_1.fastq.gz | 29MB | 272,479 | R1 (forward) |
| sample2_1_trimmed_2.fastq.gz | 32MB | 272,479 | R2 (reverse) |
| chr22.fasta | 49MB | 1 sequence | Reference genome |
| chr22.fasta.* | 255MB | N/A | BWA-MEM2 index |

**Total data size:** ~365MB (includes indexes)

## Performance Results Summary

### Method Comparison (on 4-core system)

| Method | Total Time | Read I/O | Alignments | Status |
|--------|-----------|----------|------------|---------|
| Internal Reader | 89.70s | 4.94s | 545,788 | ✅ Baseline |
| **Pigz Streaming** | **82.58s** | **2.28s** | **545,788** | 🏆 **Winner** |
| Fixed FASTR Streaming | 95.43s | 13.8s* | 545,788 | ✅ Works |

*Includes conversion overhead

### Key Findings

1. **Pigz is 7.9% faster** than internal reader
   - 54% faster decompression (2.28s vs 4.94s)
   - Better CPU utilization (4 BWA + 4 pigz threads)

2. **Fixed FASTR works correctly** 
   - All 272,479 reads preserved
   - Proper paired-end matching
   - Identical alignment results
   - But slower due to conversion overhead

3. **FASTR best for archival**
   - 21% better compression than FASTQ.GZ
   - Not suitable for active analysis pipelines

## Running Your Own Benchmarks

### Modify Thread Count

Edit the scripts to test different thread configurations:

```bash
# simple_benchmark.sh, line 22:
THREADS=8  # Change from 4 to 8

# Adjust pigz threads (line 58-59):
# -p 4  # Change from 2 to 4 per file (total 8 threads)
```

### Use Different Test Data

Replace the test data files:

```bash
cp your_reads_1.fastq.gz app/bwa-mem2/test_data/
cp your_reads_2.fastq.gz app/bwa-mem2/test_data/
cp your_reference.fasta app/bwa-mem2/test_data/

# Rebuild BWA index if needed:
cd app/bwa-mem2/test_data/
bwa-mem2 index your_reference.fasta
```

### Profile Performance

To get detailed timing information:

```bash
# Add -v flag for verbose output
pixi run benchmark-simple
# Output includes detailed timing breakdowns
```

## Troubleshooting

### Issue: "pigz: command not found"

**Solution:** Ensure pixi environment is activated
```bash
pixi install  # Reinstall dependencies
pixi run benchmark-simple  # Use pixi run
```

### Issue: "python: command not found"

**Solution:** Activate the pixi environment
```bash
pixi shell
# Now you're in the environment
python --version
bash performance/simple_benchmark.sh
```

### Issue: Test data files missing

**Solution:** Check if files exist and permissions are correct
```bash
ls -lh app/bwa-mem2/test_data/
# If empty, verify they were copied:
cd FASTR/app/bwa-mem2
cp ../../test_data/* test_data/  # Adjust path as needed
```

### Issue: BWA-MEM2 index not found

**Solution:** The index should be included. Rebuild if needed:
```bash
cd app/bwa-mem2/test_data/
bwa-mem2 index chr22.fasta
cd ..
```

### Issue: Alignment produces 0 alignments

**Solution:** Check if FASTR reconstruction was successful
```bash
# Verify read counts match between R1 and R2
pixi run verify-pairing

# Check FASTR file integrity
python ../../src/toFASTQ.py app/bwa-mem2/test_data/sample2_1_trimmed_1_fixed.fastr /tmp/test.fastq --threads 1
samtools view -c test.bam  # Compare alignment counts
```

## Development Notes

### For FASTR Developers

If modifying chunk processing:

1. **Always test with `verify-pairing` first**
   - Ensures input data is valid
   - Confirms pairing synchronization

2. **Run `benchmark-fastr` after changes**
   - Tests the full pipeline
   - Verifies no data loss
   - Checks alignment accuracy

3. **Compare with baseline methods**
   - Run `benchmark-simple` 
   - Ensure comparable or better performance

### For Contributing Bug Fixes

1. **Create a test case** that reproduces the bug
2. **Add a script** in `app/bwa-mem2/` directory
3. **Document the issue** in comments
4. **Verify the fix** with all benchmarks
5. **Add timing comparison** before/after

## Environment Specifications

### Tested On

- **OS:** Linux (Ubuntu 20.04+)
- **Python:** 3.9 - 3.12
- **CPU Arch:** x86_64 (AVX2 support)
- **RAM:** 8GB minimum (16GB recommended)

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| python | >=3.9 | FASTR implementation |
| numpy | >=1.24 | Array operations |
| numba | >=0.57 | JIT compilation for speed |
| bwa-mem2 | >=2.2.1 | Read alignment |
| samtools | >=1.18 | BAM file processing |
| pigz | >=2.7 | Parallel gzip decompression |

## References

- **FASTR Repository:** https://github.com/ALSER-Lab/FASTR
- **BWA-MEM2:** https://github.com/bwa-mem2/bwa-mem2
- **Pixi:** https://pixi.sh

## Continuous Integration / GitHub Actions

This project includes automated testing via GitHub Actions. The workflow verifies the FASTR bug fix and data integrity on every push and pull request.

### Workflow: `test.yml`

**Triggers:** Push to `main` or `develop` branches, Pull requests to `main` or `develop`

**Steps:**
1. **Set up Pixi** - Installs Pixi environment manager (with caching)
2. **Cache test data** - Caches test files in `app/bwa-mem2/test_data/` to speed up repeated runs
3. **Install dependencies** - Sets up Python 3.11 with all required packages
4. **Verify FASTQ pairing** - Confirms input files are properly synchronized
5. **Test FASTR conversion** - Runs the fixed chunk processor to verify all reads are preserved
6. **Run unit tests** - Attempts to run pytest (optional if no tests exist)

### Test Data Caching

The workflow uses GitHub's action cache to store test data in `app/bwa-mem2/test_data/` between runs.

**Cache key:** `fastr-test-data-${{ hashFiles('app/bwa-mem2/test_data/*') }}`

### Running Workflows Locally

To test the workflow locally before pushing to GitHub, use `act`:

```bash
# List available workflows
act --list

# Run the test workflow
act push -j test

# Run specific step only
act push -j test --step "Verify FASTQ pairing"
```

### Debugging CI/CD Failures

If the workflow fails:

1. Check the detailed logs in GitHub Actions tab
2. Run locally with `act` to reproduce: `act push -j test`
3. Test individual steps: `pixi run verify-pairing`, `pixi run verify-fastr-conversion`
4. Check that test data files exist: `ls -lh app/bwa-mem2/test_data/`

## Contributing

To contribute improvements to performance testing:

1. Fork the FASTR repository
2. Create a branch: `git checkout -b performance/your-improvement`
3. Add or modify benchmark scripts in `app/bwa-mem2/`
4. Update this README with new tasks
5. Test thoroughly: `pixi run benchmark-all`
6. Submit a pull request with results

## License

This performance testing setup is part of the FASTR project and follows the same license.

---

**Last Updated:** 2026-03-18  
**FASTR Version:** Main branch (with chunk processor fix)  
**Status:** ✅ Ready for production use
