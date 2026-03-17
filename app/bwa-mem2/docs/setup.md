# FASTR Performance Testing Setup - Quick Reference

## Setup Checklist

Use this guide to verify your FASTR performance testing environment is correctly configured.

### ✅ Pre-Installation

- [ ] Running Linux or macOS with bash
- [ ] Git installed: `git --version`
- [ ] ~500MB free disk space (for test data and results)
- [ ] 8GB+ RAM for running benchmarks

### ✅ Installation Steps

```bash
# 1. Install Pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# 2. Clone FASTR
git clone https://github.com/ALSER-Lab/FASTR.git
cd FASTR

# 3. Install environment
pixi install

# 4. Verify installation
pixi run verify-pairing
```

### ✅ File Structure Verification

After setup, your FASTR directory should look like:

```
FASTR/
├── pixi.toml                    # ✅ Pixi environment config
├── README.developer.md          # ✅ This guide
├── performance/
│   ├── test_data/
│   │   ├── sample2_1_trimmed_1.fastq.gz    # ✅ Test data (29MB)
│   │   ├── sample2_1_trimmed_2.fastq.gz    # ✅ Test data (32MB)
│   │   ├── chr22.fasta                     # ✅ Reference (49MB)
│   │   └── chr22.fasta.*                   # ✅ Index files (255MB)
│   ├── results/                            # Output directory
│   ├── verify_pairing.sh                   # ✅ Validation script
│   ├── simple_benchmark.sh                 # ✅ Benchmark script
│   ├── test_fastr_fixed.sh                 # ✅ FASTR test script
│   ├── BENCHMARK_RESULTS.md                # ✅ Results doc
│   ├── FASTR_BUG_FIX_SUMMARY.md           # ✅ Bug fix doc
│   └── FINAL_RESULTS_WITH_FIXED_FASTR.md  # ✅ Final results doc
├── src/
│   ├── toFASTR.py                          # FASTQ → FASTR
│   ├── toFASTQ.py                          # FASTR → FASTQ
│   └── toFASTR_chunk_processor.py          # ✅ FIXED
└── ... other files ...
```

### ✅ Dependency Check

Verify all dependencies are installed:

```bash
pixi run python --version          # Python 3.9+
pixi run python -c "import numpy"  # NumPy
pixi run python -c "import numba"  # Numba
pixi run bwa-mem2 -h               # BWA-MEM2
pixi run samtools --version        # SAMtools
pixi run pigz --version            # Pigz
```

### ✅ Test Data Validation

```bash
cd FASTR/performance

# Check file sizes
ls -lh test_data/

# Verify FASTQ files are readable
zcat test_data/sample2_1_trimmed_1.fastq.gz | head -4
zcat test_data/sample2_1_trimmed_2.fastq.gz | head -4

# Verify BWA index
ls -lh test_data/chr22.fasta.*
```

### ✅ Quick Test Run

```bash
cd FASTR

# 1. Verify pairing (1-2 minutes)
pixi run verify-pairing

# Expected output:
# R1 reads: 272479
# R2 reads: 272479
# ✅ Read counts match
```

---

## Common Setup Issues & Solutions

### Issue 1: "pixi: command not found"

**Problem:** Pixi not in PATH

**Solution:**
```bash
# Add to ~/.bashrc or ~/.zshrc:
export PATH="$HOME/.pixi/bin:$PATH"

# Reload shell:
source ~/.bashrc  # or source ~/.zshrc
```

### Issue 2: "Disk space: 120GB available, but only 50MB test data"

**Problem:** Test data files are missing

**Solution:**
```bash
# Check what's in performance/test_data/
ls -lh FASTR/performance/test_data/

# If empty, files need to be copied from original location
# Example (adjust path as needed):
cp /path/to/original/test/data/* FASTR/performance/test_data/
```

### Issue 3: "python: command not found" or version mismatch

**Problem:** Python from pixi environment not found

**Solution:**
```bash
# Always use pixi to run Python commands:
pixi run python --version

# Or activate the environment:
pixi shell
python --version
exit
```

### Issue 4: "bwa-mem2: command not found"

**Problem:** BWA-MEM2 not installed

**Solution:**
```bash
# Reinstall environment:
cd FASTR
pixi install

# Verify:
pixi run bwa-mem2 -h
```

### Issue 5: "Test produced 0 alignments"

**Problem:** Could indicate FASTR data loss or file issues

**Solution:**
```bash
# First, verify input files:
pixi run verify-pairing

# If that passes, check FASTR reconstruction:
pixi run python
>>> import subprocess
>>> result = subprocess.run([
...     'python', 'src/toFASTQ.py',
...     'performance/test_data/sample2_1_trimmed_1.fastq',
...     '/tmp/test.fastq', '--threads', '1'
... ], capture_output=True)
>>> print(result.returncode)  # Should be 0
```

---

## Running Benchmarks

### Minimal Test (1-2 minutes)

```bash
pixi run verify-pairing
```

### Full Benchmark Suite (~8 minutes)

```bash
# Option 1: Run all at once
pixi run benchmark-all

# Option 2: Run individually
pixi run benchmark-simple        # 3 minutes
pixi run benchmark-fastr         # 5 minutes
```

### Individual Benchmarks

```bash
# Compare internal reader vs pigz
pixi run benchmark-simple

# Test FASTR format (fixed version)
pixi run benchmark-fastr

# Validate input files
pixi run verify-pairing
```

---

## Understanding the Bug Fix

### The Problem (Before Fix)

```
FASTQ → FASTR Conversion
Input: 272,479 reads
↓ (buggy code drops reads at chunk boundaries)
Output: 272,472 reads  ❌ 7 READS LOST!
```

**Root cause:** Chunk boundary detection using `buffer.rfind(b"\n@")` fails because quality scores can contain "@" characters.

### The Solution (After Fix)

```
FASTQ → FASTR Conversion
Input: 272,479 reads
↓ (validates record structure properly)
Output: 272,479 reads  ✅ ALL READS PRESERVED
```

**Fix:** Proper FASTQ record validation that checks:
- Header starts with "@"
- Separator starts with "+"
- Quality length matches sequence length

### Verification

Run this to confirm the fix works:

```bash
pixi run benchmark-fastr

# Look for this in output:
# Total sequences reconstructed: 272,479
# 545788 + 0 in total (QC-passed reads)
```

---

## Performance Results

### Quick Summary

| Method | Time | Speed | Status |
|--------|------|-------|--------|
| Internal Reader | 89.70s | Baseline | ✅ |
| **Pigz Streaming** | **82.58s** | **+7.9%** | 🏆 |
| FASTR (Fixed) | 95.43s | -6.4% | ✅ |

### What This Means

- **Use pigz streaming** for active analysis pipelines (fastest)
- **Use FASTR** for archival storage (21% better compression)
- **All three methods** now produce identical alignment results

---

## Next Steps

### After Successful Setup

1. **Review the bug fix:**
   ```bash
   # Read the fix summary
   cat FASTR/performance/FASTR_BUG_FIX_SUMMARY.md
   
   # Check the fixed code
   nano FASTR/src/toFASTR_chunk_processor.py
   # Look for: find_last_fastq_record_boundary()
   ```

2. **Run full benchmarks:**
   ```bash
   pixi run benchmark-all
   
   # View results
   cd FASTR/performance/results/
   samtools flagstat test_fastr_fixed_streaming.bam
   ```

3. **Test with your own data:**
   ```bash
   # Copy your FASTQ files to performance/test_data/
   # Update the scripts to use your files
   # Run pixi run benchmark-simple
   ```

### Contributing Improvements

Found an issue or want to improve the benchmarks?

1. Create a feature branch: `git checkout -b perf/your-improvement`
2. Modify scripts in `FASTR/performance/`
3. Test thoroughly: `pixi run benchmark-all`
4. Document changes in `README.developer.md`
5. Submit a pull request

---

## Reference Documents

Located in `FASTR/performance/`:

- **`FASTR_BUG_FIX_SUMMARY.md`** - Detailed technical analysis of the bug and fix
- **`FINAL_RESULTS_WITH_FIXED_FASTR.md`** - Complete benchmark results with all three methods
- **`BENCHMARK_RESULTS.md`** - Initial benchmarks (pre-fix)

---

## Support

### Debugging

If benchmarks fail, check:

```bash
# 1. Pixi environment
pixi info

# 2. Available tools
pixi run which bwa-mem2
pixi run which samtools
pixi run which pigz

# 3. Test data
ls -lh FASTR/performance/test_data/ | wc -l  # Should be 9 files

# 4. FASTR installation
pixi run python -c "import sys; sys.path.insert(0, 'src'); from toFASTR_chunk_processor import find_last_fastq_record_boundary; print('OK')"
```

### Getting Help

1. Check this guide's "Common Setup Issues" section
2. Review `README.developer.md` for detailed information
3. Check log files in `FASTR/performance/results/`
4. Open an issue at: https://github.com/ALSER-Lab/FASTR/issues

---

**Status:** ✅ Ready to use  
**Last Updated:** 2026-03-18  
**Version:** FASTR with fixed chunk processor
