# FASTQ Compression Benchmark
Reproducible scripts for benchmarking compression and decompression of FASTQ files across general-purpose, fastq-specific  and fastq-alternative tools.

Compatible with any FASTQ-like text input.

---

## Overview

These scripts compresses input FASTQ, decompresses the file and logs the run details to disk which is also appended as metrics to a TSV file.

All measurements use `/usr/bin/time` for wall-clock timing and peak memory (Max RSS). 
```
.
├── scripts/
│   ├── spring.sh
│   ├── gzip.sh
│   ├── bzip2.sh
│   ├── xz.sh
│   ├── bsc.sh
│   ├── zstd.sh
│   ├── zip.sh
│   ├── pigz.sh
│   ├── renano.sh
│   └── sam.sh
├── outputs/
│   └── <tool>-output/
│       ├── logs/
│       │   ├── metrics.tsv
│       │   └── <tool>_<timestamp>.log
│       └── <run artifacts>
└── README.md
```
---

## Metrics Schema

Each run appends one row to `metrics.tsv`:
| Column | Description |
|--------|-------------|
| `tool` | Tool identifier (e.g., `spring`, `gzip`) |
| `run_id` | Timestamp (`YYYYMMDD_HHMMSS`) |
| `threads_requested` | Threads requested (even if tool ignores it) |
| `thread_flag_compress` | Flag used for compression threads (`-t`, `-T`, `-@`, or `NA`) |
| `thread_flag_decompress` | Flag used for decompression threads |
| `mt_compress` | `1` if compression supports multithreading, else `0` |
| `mt_decompress` | `1` if decompression supports multithreading, else `0` |
| `compression_time_s` | Wall time for compression (seconds) |
| `compressed_size_mb` | Compressed output size (MiB) |
| `compression_ratio` | `input_raw_bytes / compressed_bytes` |
| `peak_mem_compress_mb` | Peak RSS during compression (MiB) |
| `decompression_time_s` | Wall time for decompression (seconds) |
| `decompressed_size_mb` | Decompressed output size (MiB) |
| `peak_mem_decompress_mb` | Peak RSS during decompression (MiB) |
| `input_matches_decompressed` | `1` if byte-stream comparison passes, else `0` |

**Note on validation:** We use `cmp -s` for strict byte-stream comparison. Tools that reorder reads (e.g., reference-based CRAM) or apply lossy compression will return `0` even if scientifically valid.

---

## Installation

### Via Conda (Recommended)
```bash
conda create -n fastq_bench -c conda-forge -c bioconda \
  spring gzip bzip2 xz zstd pigz bsc zip unzip samtools picard

conda activate fastq_bench
```