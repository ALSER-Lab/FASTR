#!/usr/bin/env bash
set -euo pipefail

# Benchmark wrapper for FASTQ -> SAM/BAM/CRAM workflows.
# Handles four main tracks:
#   1. bam            (FastqToSam -> unaligned BAM)
#   2. sam            (FastqToSam -> BAM -> SAM)
#   3. cram_unaligned (FastqToSam -> BAM -> CRAM)
#   4. cram_aligned   (minimap2 -> CRAM)
#
# Usage:
#   ./samtools.sh --input1 reads.fastq --format bam --threads 16 --tmpdir /tmp --out-prefix /path/out/HG002_run1_bam
#   ./samtools.sh --input1 reads.fastq --format sam --threads 16 --tmpdir /tmp --out-prefix /path/out/HG002_run2_sam
#   ./samtools.sh --input1 reads.fastq --format cram_unaligned --threads 16 --tmpdir /tmp --out-prefix /path/out/HG002_run3_cram_unaligned
#   ./samtools.sh --input1 reads.fastq --format cram_aligned --threads 16 --tmpdir /tmp --out-prefix /path/out/HG002_run4_cram_aligned --ref /path/ref.fa
#

#
# Defaults
#
IN1=""
IN2=""
FORMAT=""
THREADS="4"
TMPDIR="/tmp"
OUT_PREFIX="samtools_output"
REF=""
SAMPLE=""

TOOL="samtools"

#
# Parse named args
#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input|--input1) IN1="$2"; shift 2;;
    --input2) IN2="$2"; shift 2;;
    --format) FORMAT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    -h|--help)
      cat <<EOF
Usage:
  $0 --input1 <fastq[.gz]> [--input2 <fastq[.gz]>] --format <bam|sam|cram_unaligned|cram_aligned> \\
     [--threads N] [--tmpdir DIR] [--out-prefix PREFIX] [--ref REF.fa] [--sample NAME]

Notes:
  - bam/sam/cram_unaligned use Picard FastqToSam (unaligned BAM) as the FASTQ container step.
  - cram_aligned uses minimap2 alignment + samtools CRAM encoding (not byte-identical FASTQ roundtrip).
EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1;;
  esac
done

[[ -n "$IN1" ]] || { echo "ERROR: --input1 is required"; exit 1; }
[[ -n "$FORMAT" ]] || { echo "ERROR: --format is required"; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }

#
# Sanity checks
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH (activate conda env)."; exit 1; }

# Picard needed for unaligned track - check
if [[ "$FORMAT" == "bam" || "$FORMAT" == "sam" || "$FORMAT" == "cram_unaligned" ]]; then
  command -v picard >/dev/null 2>&1 || { echo "ERROR: picard not found in PATH (needed for bam/sam/cram_unaligned). Install/activate picard from conda."; exit 1; }
fi

# minimap2 needed for aligned track - check
if [[ "$FORMAT" == "cram_aligned" ]]; then
  command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not found in PATH (needed for cram_aligned)."; exit 1; }
  [[ -n "$REF" ]] || { echo "ERROR: --ref is required for --format cram_aligned"; exit 1; }
  [[ -f "$REF" ]] || { echo "ERROR: Reference not found: $REF"; exit 1; }
  [[ -f "${REF}.fai" ]] || echo "WARNING: ${REF}.fai not found. Indexing may be attempted (can be slow / can fail if ref dir not writable)."
fi

#
# Directories checks
#
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

OUT_PREFIX="${OUT_PREFIX%/}"
OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"

# Logging
LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${FORMAT}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"
# Pipe everything to log
exec > >(tee -a "$LOGFILE") 2>&1

# Sample default
if [[ -z "$SAMPLE" ]]; then
  SAMPLE="$BASE_NAME"
fi

echo "=== samtools benchmark ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "format=$FORMAT"
echo "run_id=$RUN_ID"
echo "input1=$IN1"
echo "input2=${IN2:-NA}"
echo "threads_requested=$THREADS"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "ref=${REF:-NA}"
echo "sample=$SAMPLE"
echo "logfile=$LOGFILE"
echo

# Helper for /usr/bin/time to capture RSS memory
run_tm() {
  local label="$1"; shift
  local tf
  tf="$(mktemp -p "$TMPDIR" "time_${label}_XXXX.txt")"

  { echo; echo "---- $label ----"; echo "CMD: $*"; echo "Start($label): $(date -Is)"; } >&2

  set +e
  /usr/bin/time -f "%e\t%M" -o "$tf" "$@" 1>&2
  local rc=$?
  set -e

  echo "End($label): $(date -Is) (rc=$rc)" >&2

  if [[ $rc -ne 0 ]]; then
    rm -f "$tf"
    return $rc
  fi
# Return seconds and MB
  awk -F'\t' 'END{printf "%.6f\t%.3f\n",$1,($2/1024.0)}' "$tf"
  rm -f "$tf"
}

# Helper to handle .gz or plain inputs
stream() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -cd -- "$f"
  else
    cat -- "$f"
  fi
}

#
# Calculate Raw input bytes (uncompressed)
#
echo "Computing raw input bytes (uncompressed)..."
if [[ -n "$IN2" ]]; then
  in1_bytes="$(stream "$IN1" | wc -c)"
  in2_bytes="$(stream "$IN2" | wc -c)"
  in_bytes="$((in1_bytes + in2_bytes))"
else
  in_bytes="$(stream "$IN1" | wc -c)"
fi
echo "input_raw_bytes=$in_bytes"
echo

# Temp workspace for intermediate files (e.g., unaligned BAM before CRAM)
WORKDIR="$OUT_DIR/.work_${BASE_NAME}_${RUN_ID}"
mkdir -p "$WORKDIR"

# Format config + unified TSV meta
MT_COMPRESS=1
MT_DECOMPRESS=1
THREAD_FLAG_DECOMPRESS="-@"  # samtools fastq always uses -@

ARCHIVE=""
DEC1=""
DEC2=""

case "$FORMAT" in
  bam)
    # Compression is Picard only (generally single-threaded)
    THREAD_FLAG_COMPRESS="picard:NA"
    MT_COMPRESS=0
    ARCHIVE="${OUT_PREFIX}.bam"
    DEC1="${OUT_PREFIX}.dec_R1.fastq"
    DEC2="${OUT_PREFIX}.dec_R2.fastq"
    ;;
  sam)
    # Picard (ST) + samtools view (MT)
    THREAD_FLAG_COMPRESS="picard:NA;samtools:-@"
    ARCHIVE="${OUT_PREFIX}.sam"
    DEC1="${OUT_PREFIX}.dec_R1.fastq"
    DEC2="${OUT_PREFIX}.dec_R2.fastq"
    ;;
  cram_unaligned)
    # Picard (ST) + samtools view -C (MT)
    THREAD_FLAG_COMPRESS="picard:NA;samtools:-@"
    ARCHIVE="${OUT_PREFIX}.cram"
    DEC1="${OUT_PREFIX}.dec_R1.fastq"
    DEC2="${OUT_PREFIX}.dec_R2.fastq"
    ;;
  cram_aligned)
    # minimap2 + samtools view (both MT)
    THREAD_FLAG_COMPRESS="minimap2:-t;samtools:-@"
    ARCHIVE="${OUT_PREFIX}.cram"
    DEC1="${OUT_PREFIX}.dec_R1.fastq"
    DEC2="${OUT_PREFIX}.dec_R2.fastq"
    ;;
  *)
    echo "ERROR: Invalid --format '$FORMAT'. Choose: bam, sam, cram_unaligned, cram_aligned"
    exit 1
    ;;
esac

echo "archive=$ARCHIVE"
echo "decompressed_R1=$DEC1"
echo "decompressed_R2=${IN2:+$DEC2}"
echo

# Compression
rm -f "$ARCHIVE"

if [[ "$FORMAT" == "cram_aligned" ]]; then
# Pipeline: minimap2 -> sort/view -> CRAM
  # Using bash -c to keep the pipe contained for timing
  read c_time c_peak_mb < <(
    run_tm "compress" bash -c '
      ref="$1"; in1="$2"; in2="$3"; out="$4"; th="$5"
      if [[ -n "$in2" ]]; then
        # Paired-end: interleave for minimap2
        minimap2 -a -x sr -t "$th" "$ref" "$in1" "$in2" \
          | samtools view -@ "$th" -C -T "$ref" -o "$out" -
      else
        minimap2 -a -x sr -t "$th" "$ref" "$in1" \
          | samtools view -@ "$th" -C -T "$ref" -o "$out" -
      fi
    ' cram_aligned_pipe "$REF" "$IN1" "${IN2:-}" "$ARCHIVE" "$THREADS"
  )
else
  # Unaligned track uses Picard FastqToSam to make BAM first
  U_BAM="$WORKDIR/unaligned.bam"
  rm -f "$U_BAM"

  # Picard FastqToSam (single-end or paired-end)
  if [[ -n "$IN2" ]]; then
    # Note: Picard FastqToSam doesn't reliably use threads; we still record threads_requested for consistency.
    read p_time p_peak_mb < <(
      run_tm "picard_fastqtosam" picard FastqToSam \
        F1="$IN1" F2="$IN2" \
        O="$U_BAM" \
        SM="$SAMPLE" \
        RG="rg1"
    )
  else
    read p_time p_peak_mb < <(
      run_tm "picard_fastqtosam" picard FastqToSam \
        F1="$IN1" \
        O="$U_BAM" \
        SM="$SAMPLE" \
        RG="rg1"
    )
  fi
  [[ -s "$U_BAM" ]] || { echo "ERROR: Picard did not produce unaligned BAM: $U_BAM"; exit 1; }

  # Now produce requested archive
  if [[ "$FORMAT" == "bam" ]]; then
    # Archive IS the BAM from Picard
    mv -f "$U_BAM" "$ARCHIVE"
    # For unified timing/memory, treat compression as picard step
    c_time="$p_time"
    c_peak_mb="$p_peak_mb"
  elif [[ "$FORMAT" == "sam" ]]; then
    # BAM -> SAM
    read s_time s_peak_mb < <(
      run_tm "bam_to_sam" samtools view -@ "$THREADS" -h -O SAM -o "$ARCHIVE" "$U_BAM"
    )
    # Compression time includes Picard + conversion (pipeline)
    c_time="$(awk -v a="$p_time" -v b="$s_time" 'BEGIN{printf "%.6f", a+b}')"
    c_peak_mb="$(awk -v a="$p_peak_mb" -v b="$s_peak_mb" 'BEGIN{printf "%.3f", (a>b)?a:b}')"
    rm -f "$U_BAM"
  elif [[ "$FORMAT" == "cram_unaligned" ]]; then
    # BAM -> unaligned CRAM (no reference)
    read s_time s_peak_mb < <(
      run_tm "bam_to_cram" samtools view -@ "$THREADS" -C -o "$ARCHIVE" "$U_BAM"
    )
    c_time="$(awk -v a="$p_time" -v b="$s_time" 'BEGIN{printf "%.6f", a+b}')"
    c_peak_mb="$(awk -v a="$p_peak_mb" -v b="$s_peak_mb" 'BEGIN{printf "%.3f", (a>b)?a:b}')"
    rm -f "$U_BAM"
  else
    echo "ERROR: Unexpected FORMAT in unaligned track: $FORMAT"
    exit 1
  fi
fi

[[ -s "$ARCHIVE" ]] || { echo "ERROR: compression did not produce archive: $ARCHIVE"; exit 1; }

# Stats
out_bytes="$(stat -c%s "$ARCHIVE")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"

# Decompression (archive -> FASTQ)
rm -f "$DEC1" "${IN2:+$DEC2}"

if [[ -n "$IN2" ]]; then
  # Paired-end output
  DECOMP_ARGS=(samtools fastq -@ "$THREADS" -n -1 "$DEC1" -2 "$DEC2" -0 /dev/null -s /dev/null)
else
  # Single-end output (R1 only)
  DECOMP_ARGS=(samtools fastq -@ "$THREADS" -n -0 "$DEC1")
fi

if [[ "$FORMAT" == "cram_aligned" ]]; then
  # reference needed to avoid fetch attempts / ensure decode
  DECOMP_ARGS+=(--reference "$REF")
fi

DECOMP_ARGS+=("$ARCHIVE")

read d_time d_peak_mb < <(run_tm "decompress" "${DECOMP_ARGS[@]}")

# Validate outputs exist
[[ -s "$DEC1" ]] || { echo "ERROR: decompression did not produce $DEC1"; exit 1; }
if [[ -n "$IN2" ]]; then
  [[ -s "$DEC2" ]] || { echo "ERROR: decompression did not produce $DEC2"; exit 1; }
fi

# Decompressed size (sum if paired)
dec_bytes="$(stat -c%s "$DEC1")"
if [[ -n "$IN2" ]]; then
  dec_bytes2="$(stat -c%s "$DEC2")"
  dec_bytes="$((dec_bytes + dec_bytes2))"
fi
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"

# Verify integrity: compare decompressed FASTQ to original input FASTQ
# Note: Byte matching almost always fails for CRAM/BAM roundtrips because headers,
# sort order, and compression blocks change. We flag it, but don't fail the run.
input_matches_decompressed=0
if [[ "$FORMAT" == "cram_aligned" ]]; then
  echo "input_matches_decompressed=0 (skipped: aligned CRAM reorders / may change FASTQ roundtrip)"
else
  if [[ -n "$IN2" ]]; then
    if cmp -s <(stream "$IN1") "$DEC1" && cmp -s <(stream "$IN2") "$DEC2"; then
      input_matches_decompressed=1
    fi
  else
    if cmp -s <(stream "$IN1") "$DEC1"; then
      input_matches_decompressed=1
    fi
  fi
  echo "input_matches_decompressed=$input_matches_decompressed"
  echo "NOTE: For HTS container conversions, FASTQ headers/order may change; cmp can be 0 even if bases/quals are preserved."
fi

# Metrics TSV 
if [[ ! -f "$METRICS_TSV" ]]; then
  printf "tool\trun_id\tthreads_requested\tthread_flag_compress\tthread_flag_decompress\tmt_compress\tmt_decompress\tcompression_time_s\tcompressed_size_mb\tcompression_ratio\tpeak_mem_compress_mb\tdecompression_time_s\tdecompressed_size_mb\tpeak_mem_decompress_mb\tinput_matches_decompressed\n" \
    >> "$METRICS_TSV"
fi

printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.6f\t%s\t%.6f\t%.3f\t%.6f\t%s\t%.3f\t%d\n" \
  "${TOOL}_${FORMAT}" "$RUN_ID" "$THREADS" "$THREAD_FLAG_COMPRESS" "$THREAD_FLAG_DECOMPRESS" \
  "$MT_COMPRESS" "$MT_DECOMPRESS" \
  "$c_time" "$comp_size_mb" "$comp_ratio" "$c_peak_mb" \
  "$d_time" "$dec_size_mb" "$d_peak_mb" \
  "$input_matches_decompressed" \
  | tee -a "$METRICS_TSV"

echo
echo "Done: $(date -Is)"
echo "Log: $LOGFILE"
echo "Metrics: $METRICS_TSV"