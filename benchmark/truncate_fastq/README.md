# Fast FASTQ File Truncation on HPC

A quick guide for truncating large FASTQ files using `dd` on a SLURM-based HPC cluster.

## Overview

This workflow uses `dd` to efficiently extract the first 33GB from a large HiFi FASTQ file. Using `dd` with a large block size is significantly faster than `head -c` for large files.

## Prerequisites

- Access to a SLURM-based HPC cluster
- Read access to the source FASTQ file
- Write access to the output directory

## Usage

```bash
sbatch -p qDEV -N 1 -n 48 -t 11:00:00 -o wordcount.txt --wrap="dd if=/path/to/input.fastq of=/path/to/output.fastq bs=1M iflag=count_bytes count=33000000000"
```

## SLURM Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-p` | qDEV | Partition/queue name |
| `-N` | 1 | Number of nodes |
| `-n` | 48 | Number of CPU cores |
| `-t` | 11:00:00 | Time limit (11 hours) |
| `-o` | wordcount.txt | Output log file |

## dd Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `if` | input.fastq | Input file path |
| `of` | output.fastq | Output file path |
| `bs` | 1M | Block size (1 megabyte) |
| `iflag=count_bytes` | - | Interpret `count` as bytes |
| `count` | 33000000000 | Number of bytes to copy (~33GB) |

## Why dd over head?

| Method | Buffer Size | Speed |
|--------|-------------|-------|
| `head -c` | 8KB (default) | Slower |
| `dd bs=1M` | 1MB | Much faster |

For large files, `dd` with a larger block size reduces the number of read/write operations significantly.

## Example

```bash
sbatch -p qDEV -N 1 -n 48 -t 11:00:00 -o truncate.log --wrap="dd if=/home/users/muddin21/WA/projects/fastr/fastq/HiFi_fastq/ERR16043555.fastq of=/home/users/muddin21/WA/projects/fastr/fastq/HiFi_fastq/ERR16043555_trunc.fastq bs=1M iflag=count_bytes count=33000000000"
```

## Monitor Job

```bash
# Check job status
squeue -u $USER

# View output log
tail -f wordcount.txt
```
