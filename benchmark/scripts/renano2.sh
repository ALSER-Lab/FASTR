#!/usr/bin/env bash
set -euo pipefail


# Usage:
#   ./renano.sh \
#     --input ./fastq_files/D1_S1_L001_R1_001-017.fastq \
#     --threads 32 \
#     --tmpdir /tmp \
#     --out-prefix /folder/HG002_R1_001-017

# Configs
IN=""
THREADS="16"
TMPDIR="/tmp"
OUT_PREFIX="renano-output"
TOOL="renano"
# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file.fastq> [--threads N] [--tmpdir DIR] [--out-prefix PREFIX]

Notes:
  - This script expects an UNCOMPRESSED .fastq input.
  - RENANO supports threads with -t for both compression and decompression.
EOF
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

[[ -n "$IN" ]] || { echo "ERROR: --input is required"; exit 1; }

# Guard against passing a directory or trailing slash
OUT_PREFIX="${OUT_PREFIX%/}"

# Sanity Checks 
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

command -v renano >/dev/null 2>&1 || { echo "ERROR: renano not found in PATH. Activate your conda env."; exit 1; }

# RENANO (per README) does not advertise .gz input; enforce uncompressed for reproducibility
if [[ "$IN" == *.gz ]]; then
  echo "ERROR: RENANO wrapper expects uncompressed .fastq (got .gz)."
  echo "       Please gunzip outside the benchmark and re-run."
  exit 1
fi

# Threading metadata (RENANO supports -t in both phases)
THREAD_FLAG_COMPRESS="-t"
THREAD_FLAG_DECOMPRESS="-t"
MT_COMPRESS=1
MT_DECOMPRESS=1

# Output paths
OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"

LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"

# Logging: redirect stdout/stderr to logfile and console
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== RENANO benchmark (single FASTQ) ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "threads_requested=$THREADS"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "logfile=$LOGFILE"
echo

RENANO_OUT="${OUT_PREFIX}.renano"
DEC_OUT="${OUT_PREFIX}.dec.fastq"

echo "compressed_file=$RENANO_OUT"
echo "decompressed_file=$DEC_OUT"
echo

run_tm() {
  # stdout: "<seconds>\t<peak_mb>"  (captured by read)
  # stderr: command output/progress (goes into log)
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

  awk -F'\t' 'END{printf "%.6f\t%.3f\n",$1,($2/1024.0)}' "$tf"
  rm -f "$tf"
}
# Compression (no reference)
rm -f "$RENANO_OUT"
read c_time c_peak_mb < <(run_tm "compress" renano -t "$THREADS" "$IN" "$RENANO_OUT")
[[ -s "$RENANO_OUT" ]] || { echo "ERROR: compression did not produce $RENANO_OUT"; exit 1; }

# Stats
in_bytes="$(stat -c%s "$IN")"
out_bytes="$(stat -c%s "$RENANO_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"


# Decompression
rm -f "$DEC_OUT"
read d_time d_peak_mb < <(run_tm "decompress" renano -d -t "$THREADS" "$RENANO_OUT" "$DEC_OUT")
[[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }

dec_bytes="$(stat -c%s "$DEC_OUT")"
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"


# Cerify input matches decompressed

input_matches_decompressed=0
if cmp -s -- "$IN" "$DEC_OUT"; then
  input_matches_decompressed=1
fi
echo "input_matches_decompressed=$input_matches_decompressed"


# Metrics TSV logs
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