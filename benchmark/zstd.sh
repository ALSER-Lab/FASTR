#!/usr/bin/env bash
set -euo pipefail

# zstd benchmark wrapper
# - compress + decompress
# - logs always saved
# - unified metrics TSV appended each run
#
# Examples:
#   ./zstd.sh \
#     --input ./fastq_files/D1_S1_L001_R1_001-017.fastq \
#     --threads 16 \
#     --fast 10 \
#     --tmpdir /tmp \
#     --out-prefix /folder/zstd-output/HG002_run1_fast
#
#   ./zstd.sh \
#     --input ./fastq_files/D1_S1_L001_R1_001-017.fastq \
#     --threads 16 \
#     --best \
#     --tmpdir /tmp \
#     --out-prefix /folder/zstd-output/HG002_run2_best

#
# Defaults
#
IN=""
THREADS="16"
TMPDIR="/tmp"
OUT_PREFIX="zstd-output"
TOOL="zstd"

MODE="fast"          # fast|best
FAST_LEVEL="10"      # used only in fast mode
BEST_LEVEL="19"      # -19 (variable, can be overridden)

# Threading metadata
THREAD_FLAG_COMPRESS="-T"
THREAD_FLAG_DECOMPRESS="-T"
MT_COMPRESS=1
MT_DECOMPRESS=1

# Avoid env-driven defaults
unset ZSTD_CLEVEL ZSTD_NBTHREADS

#
# Parse named args
#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;

    --fast) MODE="fast"; FAST_LEVEL="${2:-10}"; shift 2;;  
    --best) MODE="best"; shift 1;;

    --best-level) BEST_LEVEL="$2"; shift 2;;               
    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file> [--threads N] [--tmpdir DIR] [--out-prefix PREFIX] [--fast N | --best]

Notes:
  - Default is --fast 10.
  - --fast N uses: zstd --fast=N -T<threads>
  - --best uses:    zstd -<best_level> -T<threads>   (default best_level=19)
  - Encode run identifier in --out-prefix (e.g., .../run1_fast10, .../run2_best19).
EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1;;
  esac
done

[[ -n "$IN" ]] || { echo "ERROR: --input is required"; exit 1; }
OUT_PREFIX="${OUT_PREFIX%/}"

# Validate integers
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ "$FAST_LEVEL" =~ ^[0-9]+$ ]] || { echo "ERROR: --fast requires a non-negative integer"; exit 1; }
[[ "$BEST_LEVEL" =~ ^[0-9]+$ ]] || { echo "ERROR: --best-level must be an integer"; exit 1; }

# Tool checks
command -v zstd >/dev/null 2>&1 || { echo "ERROR: zstd not found in PATH (activate conda env)."; exit 1; }

# TMPDIR checks
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

# Output paths + logs
OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"

LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"

# Log
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== zstd benchmark (single file) ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "threads_requested=$THREADS"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
if [[ "$MODE" == "fast" ]]; then
  echo "mode=fast (zstd --fast=$FAST_LEVEL)"
else
  echo "mode=best (zstd -$BEST_LEVEL)"
fi
echo "logfile=$LOGFILE"
echo

ZST_OUT="${OUT_PREFIX}.zst"
DEC_OUT="${OUT_PREFIX}.dec"

echo "compressed_file=$ZST_OUT"
echo "decompressed_file=$DEC_OUT"
echo

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
  if [[ $rc -ne 0 ]]; then rm -f "$tf"; return $rc; fi

  awk -F'\t' 'END{printf "%.6f\t%.3f\n",$1,($2/1024.0)}' "$tf"
  rm -f "$tf"
}

#
# Compression
#
rm -f "$ZST_OUT"

if [[ "$MODE" == "fast" ]]; then
  read c_time c_peak_mb < <(
    run_tm "compress" bash -c 'zstd --fast="$1" -T"$2" -c -- "$3" > "$4"' \
      zstd_c "$FAST_LEVEL" "$THREADS" "$IN" "$ZST_OUT"
  )
else
  read c_time c_peak_mb < <(
    run_tm "compress" bash -c 'zstd -"$1" -T"$2" -c -- "$3" > "$4"' \
      zstd_c "$BEST_LEVEL" "$THREADS" "$IN" "$ZST_OUT"
  )
fi

[[ -s "$ZST_OUT" ]] || { echo "ERROR: compression did not produce $ZST_OUT"; exit 1; }

# Compressed size + ratio
in_bytes="$(stat -c%s "$IN")"
out_bytes="$(stat -c%s "$ZST_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"

#
# Decompression (same -T threads)
#
rm -f "$DEC_OUT"

read d_time d_peak_mb < <(
  run_tm "decompress" bash -c 'zstd -d -T"$1" -c -- "$2" > "$3"' \
    zstd_d "$THREADS" "$ZST_OUT" "$DEC_OUT"
)

[[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }

dec_bytes="$(stat -c%s "$DEC_OUT")"
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"

#
# Input matches decompressed
input_matches_decompressed=0
if cmp -s -- "$IN" "$DEC_OUT"; then
  input_matches_decompressed=1
fi
echo "input_matches_decompressed=$input_matches_decompressed"

#
# Metrics TSV (unified schema)
#
if [[ ! -f "$METRICS_TSV" ]]; then
  printf "tool\trun_id\tthreads_requested\tthread_flag_compress\tthread_flag_decompress\tmt_compress\tmt_decompress\tcompression_time_s\tcompressed_size_mb\tcompression_ratio\tpeak_mem_compress_mb\tdecompression_time_s\tdecompressed_size_mb\tpeak_mem_decompress_mb\tinput_matches_decompressed\n" \
    >> "$METRICS_TSV"
fi

printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.6f\t%s\t%.6f\t%.3f\t%.6f\t%s\t%.3f\t%d\n" \
  "$TOOL" "$RUN_ID" "$THREADS" "$THREAD_FLAG_COMPRESS" "$THREAD_FLAG_DECOMPRESS" \
  "$MT_COMPRESS" "$MT_DECOMPRESS" \
  "$c_time" "$comp_size_mb" "$comp_ratio" "$c_peak_mb" \
  "$d_time" "$dec_size_mb" "$d_peak_mb" \
  "$input_matches_decompressed" \
  | tee -a "$METRICS_TSV"

echo
echo "Done: $(date -Is)"
echo "Log saved to: $LOGFILE"
echo "Metrics appended to: $METRICS_TSV"