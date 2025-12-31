# FASTR
<p align="center">
    <a href="https://github.com/ALSER-Lab/FASTR/blob/main/figs/FASTR-Logo.png" target="_blank"><img src="https://raw.githubusercontent.com/ALSER-Lab/FASTR/refs/heads/main/figs/FASTR-Logo.png?token=GHSAT0AAAAAAC4IEXWHXJAVGQMMFZHZUGZQ2KUVFBA" alt="FASTR logo, CanvaAI-generated" width="150" border="10" /></a>
</p>

FASTR, an efficient file format designed for lossless storage of sequencing data as scalar (numerical) formats. 
FASTR transforms both textual DNA/RNA data (i.e., FASTQ) and their base quality scores into efficient/compact integer-based or binary representations.

## Features
  * FASTR is at least 2x less in size than FASTQ, and hence better to read, process, transfer.
  * FASTR can be further compressed using general-purpose compression tools, such as gzip, pigz, ...
  * Extremely fast (multithreaded) and lossless FASTR-to-FASTQ & FASTQ-to-FASTR conversion.
  * FASTR supports data from all prominent sequencing technologies (Illumina, ONT, PacBio's HiFi, and PacBio's CLR), single-end and paired-end reads, and SRA formats (https://www.ncbi.nlm.nih.gov/sra).
  * FASTR supports all widely-used Phred quality scores (Phred42, Phred63, Phred68Solexa, Phred94, Illumina RTA3, Illumina RTA4, and custom mathematical formulas).
  * Flexible Output: binar (1 uint8 per FASTR value), integer (3 uint8s per FASTR value), with/without header.
  * FASTR is compatible with minimap2 with no (or <2%) overhead, and with machine learning pipelines (i.e., numerical vectors).

<p align="center">
    <a href="https://github.com/ALSER-Lab/FASTR/blob/main/figs/FASTR-fig1.png" target="_blank"><img src="https://raw.githubusercontent.com/ALSER-Lab/FASTR/refs/heads/main/figs/FASTR-fig1.png?token=GHSAT0AAAAAAC4IEXWHEQL5XFVPABXOHOV42KUVEWQ" alt="FASTR logo, CanvaAI-generated" width="700" border="10" /></a>
</p>

-----

## Installation

Ensure you have Python 3.x installed. The tool relies on `numpy` for efficient array handling.

```bash
git clone https://github.com/ALSER-Lab/FASTR.git
cd FASTR/src
python main.py --help
```

-----

## Usage

### Basic Command

```bash
python main.py [INPUT_PATH] [OUTPUT_PATH] [FLAGS]
```

### Quick Start Modes

The `--mode` argument provides a shorthand for setting complex flag combinations.

| Mode | ID | Description | Logic Applied |
| :--- | :---: | :--- | :--- |
| **Header Only** | `1` | Header compression only. | `compress_headers=1`, `binary_write=0` |
| **Bases Only** | `2` | Base conversion into numbers only. | `compress_headers=0`, `binary_write=1` |
| **Standard** | `3` | Header + Base conversion (2 lines). | `compress_headers=1`, `binary_write=1` |
| **Compact** | `4` | Removes repeating header metadata; 1 line output. | `compress_headers=1`, `remove_repeating=1`, `binary_write=1` |

-----

## Detailed Configuration

### 1\. Paired-End & SRA

Handle paired-end sequencing data or specific SRA accessions.

| Flag | Type | Default | Description |
| :--- | :---: | :--- | :--- |
| `--paired_end` | int | 0 | Set to `1` if input contains paired-end reads. |
| `--paired_end_mode` | str | `same_file` | Output strategy: `same_file` or `separate_files`. |
| `--sra_accession` | str | None | Explicit SRA accession number (e.g., SRR12345678). |

### 2\. Header Compression

Options for reducing the footprint of sequence identifiers.

| Flag | Type | Default | Description |
| :--- | :---: | :--- | :--- |
| `--compress_headers` | int | 0 | Enable on-the-fly header compression (0/1). |
| `--sequencer_type` | str | `none` | Optimization target: `illumina`, `pacbio`, `ont`, `srr`. |
| `--multiple_flowcells` | int | 0 | Enable detection/tracking of multiple flowcells. |
| `--remove_repeating_header`| int | 0 | Strip repeating metadata, storing it only once at the top of the file. |

### 3\. Base Mapping (Grayscale)

Map nucleotides to specific integer values (0-255) for "DNA-as-Image" representations.

| Nucleotide | Flag | Default Value |
| :--- | :--- | :--- |
| **N** | `--gray_N` | 1 |
| **A** | `--gray_A` | 63 |
| **T** | `--gray_T` | 127 |
| **C** | `--gray_C` | 191 |
| **G** | `--gray_G` | 255 |

### 4\. Quality Score Processing

Manipulate Phred quality scores using standard or custom math.

| Flag | Type | Default | Description |
| :--- | :---: | :--- | :--- |
| `--extract_quality` | int | 1 | Extract quality scores (0/1). |
| `--phred_offset` | int | 33 | The Phred offset (usually 33 or 64). |
| `--min_quality` | int | 0 | Threshold for minimum quality score. |
| `--phred_alphabet` | str | `phred42` | Input alphabet: `phred42`, `phred63`, `phred94`. |
| `--quality_scaling` | str | `none` | Scaling method: `log`, `log_custom`, `custom`. |
| `--custom_formula` | str | None | **Required** if scaling is `custom`. Use `x` as the variable. |
| `--log_a` | int | 0 | Tunable variable `a` for the `log_custom` function. |

**Custom Formula Examples:**

```bash
# Linear scaling
--quality_scaling custom --custom_formula "1 + 62 * (x - 40) / 53"

# Logarithmic scaling
--quality_scaling custom --custom_formula "ln(x) * 10"
```

### 5\. Output Format

Control how the data is written to disk.

| Flag | Type | Default | Description |
| :--- | :---: | :--- | :--- |
| `--binary_write` | int | 1 | Write sequence integers as binary (1) or text (0). |
| `--keep_bases` | int | 0 | Return raw ASCII bases (A, T, C, G) without scaling. |
| `--binary_bases` | int | 0 | Use compact binary encoding for bases (A=0, T=1, etc). |
| `--keep_quality` | int | 0 | Keep original quality scores in output. |
| `--binary_quality` | int | 0 | Write quality scores as binary numeric values. |

> **Note:** You cannot use `--keep_bases 1` and `--binary_bases 1` simultaneously.

### 6\. Performance

Optimize the script for your hardware.

| Flag | Type | Default | Description |
| :--- | :---: | :--- | :--- |
| `--num_workers` | int | 1 | Number of parallel workers. **Recommended 4+ for files \>5GB.** |
| `--profile` | int | 0 | Enable `cProfile` to debug performance bottlenecks. |

-----

## Examples

**1. Quick Compression (Header + Base conversion):**

```bash
python main.py input.fastq output.bin --mode 3
```

**2. Prepare for Machine Learning (Grayscale + Custom Quality):**

```bash
python main.py data.fastq train_data.bin \
  --mode 2 \
  --quality_scaling custom \
  --custom_formula "x / 40.0" \
  --gray_A 50 --gray_T 100 --gray_C 150 --gray_G 200
```

**3. Large File Processing (Multiprocessing):**

```bash
python main.py big_genome.fastq output.bin --mode 4 --num_workers 8
```

-----
