#!/usr/bin/env bash
set -euo pipefail

# pigz benchmark wrapper 
# - compress + decompress
# - logs saved
# - unified metrics TSV appended each run
#
# Flags:
#   --fast  => pigz --fast  (equiv to -1)
#   --best  => pigz --best  (equiv to -9)
# Default: --fast
#
# Examples:
#   ./pigz.sh --input ./reads.fastq --threads 32 --fast --tmpdir /tmp \
#     --out-prefix /folder/pigz-output/run1_fast
#
#   ./pigz.sh --input ./reads.fastq --threads 32 --best --tmpdir /tmp \
#     --out-prefix /folder/pigz-output/run2_best

#
# Defaults
#
IN=""
THREADS="16"
TMPDIR="/tmp"
OUT_PREFIX="pigz-output"
TOOL="pigz"
MODE="fast"  # fast|best

# Threading metadata (what we record in TSV)
# pigz uses -p for compression; decompression is not truly parallel, but we still pass -p for consistency.
THREAD_FLAG_COMPRESS="-p"
THREAD_FLAG_DECOMPRESS="-p"
MT_COMPRESS=1
MT_DECOMPRESS=0  # per pigz docs: decompression can't be parallelized but uses helper threads
#
# Parse named args
#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --fast) MODE="fast"; shift 1;;
    --best) MODE="best"; shift 1;;
    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file.fastq[.gz]> [--threads N] [--tmpdir DIR] [--out-prefix PREFIX] [--fast|--best]

Notes:
  - Default is --fast.
  - pigz supports --fast/--best. 
  - Compression uses -p <threads>.
  - Decompression isn't truly parallelized , but pigz may use helper threads. 
  - Encode run identifier in --out-prefix (e.g., .../run1_fast, .../run2_best).
EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1;;
  esac
done

[[ -n "$IN" ]] || { echo "ERROR: --input is required"; exit 1; }
OUT_PREFIX="${OUT_PREFIX%/}"

[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

command -v pigz >/dev/null 2>&1 || { echo "ERROR: pigz not found in PATH (activate conda env)."; exit 1; }

PIGZ_FLAG="--fast"
if [[ "$MODE" == "best" ]]; then
  PIGZ_FLAG="--best"
fi

# Output paths + logs
OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"

LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"

# Log everything
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== pigz benchmark (single file) ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "threads_requested=$THREADS"
echo "mode=$MODE (flag: $PIGZ_FLAG)"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "logfile=$LOGFILE"
echo

PIGZ_OUT="${OUT_PREFIX}.gz"
DEC_OUT="${OUT_PREFIX}.dec"

echo "compressed_file=$PIGZ_OUT"
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

# Helper to stream possibly-gz input for cmp
stream() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -cd -- "$f"
  else
    cat -- "$f"
  fi
}

#
# Compression
#
rm -f "$PIGZ_OUT"
read c_time c_peak_mb < <(
  run_tm "compress" bash -c \
    'env -u GZIP -u PIGZ pigz '"$PIGZ_FLAG"' -p "$1" -c -- "$2" > "$3"' \
    pigz_c "$THREADS" "$IN" "$PIGZ_OUT"
)
[[ -s "$PIGZ_OUT" ]] || { echo "ERROR: compression did not produce $PIGZ_OUT"; exit 1; }

in_bytes="$(stat -c%s "$IN")"
out_bytes="$(stat -c%s "$PIGZ_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"

#
# Decompression (use same -p for consistency)
#
rm -f "$DEC_OUT"
read d_time d_peak_mb < <(
  run_tm "decompress" bash -c \
    'env -u GZIP -u PIGZ pigz -d -p "$1" -c -- "$2" > "$3"' \
    pigz_d "$THREADS" "$PIGZ_OUT" "$DEC_OUT"
)
[[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }

dec_bytes="$(stat -c%s "$DEC_OUT")"
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"

#
# Input matches decompressed? #
input_matches_decompressed=0
if cmp -s <(stream "$IN") <(cat -- "$DEC_OUT"); then
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