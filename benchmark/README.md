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

