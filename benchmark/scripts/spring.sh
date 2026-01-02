#!/usr/bin/env bash
set -euo pipefail

# SPRING benchmark wrapper
# Usage:
#   ./spring.sh --input file.fastq.gz --threads 16 --out-prefix ./results/run1
# Defaults
IN=""
THREADS="16"
TMPDIR="/tmp"
OUT_PREFIX="spring-output"
OUTPUT_GZIP_LEVEL="6"   # Only used if re-compressing to .gz during decompression
TOOL="spring"

# Argument parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --output-gzip-level) OUTPUT_GZIP_LEVEL="$2"; shift 2;;
    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file.fastq[.gz]> [--threads N] [--tmpdir DIR] [--out-prefix PREFIX] [--output-gzip-level N]

Notes:
  - If --input ends with .gz, SPRING is run in gzip mode (-g) and the decompressed output will be .gz.
  - Metrics TSV includes a yes/no flag: input_matches_decompressed
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

# Sanity checks
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

# Ensure SPRING exists in environment
command -v spring >/dev/null 2>&1 || { echo "ERROR: spring not found in PATH. Activate your conda env."; exit 1; }

# Threading metadata for TSV output
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

# Log everything to console + file
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== SPRING benchmark (single FASTQ) ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "threads_requested=$THREADS"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "logfile=$LOGFILE"
echo

# gzip mode based on input extension
GFLAG=""
if [[ "$IN" == *.gz ]]; then
  GFLAG="-g"
fi
echo "gzip_mode_flag=$GFLAG"

SPRING_OUT="${OUT_PREFIX}.spring"
DEC_OUT="${OUT_PREFIX}.dec.fastq"
if [[ "$GFLAG" == "-g" ]]; then
  DEC_OUT="${DEC_OUT}.gz"
fi

echo "compressed_file=$SPRING_OUT"
echo "decompressed_file=$DEC_OUT"
echo

# Timer helper (RSS memory via /usr/bin/time)
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
  # Returns: seconds \t MB
  awk -F'\t' 'END{printf "%.6f\t%.3f\n",$1,($2/1024.0)}' "$tf"
  rm -f "$tf"
}

# Compression
read c_time c_peak_mb < <(
  run_tm "compress" spring -c -i "$IN" -o "$SPRING_OUT" $GFLAG \
    -t "$THREADS" -w "$TMPDIR"
)
[[ -s "$SPRING_OUT" ]] || { echo "ERROR: compression did not produce $SPRING_OUT"; exit 1; }

# Stats
in_bytes="$(stat -c%s "$IN")"
out_bytes="$(stat -c%s "$SPRING_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"

# Decompression
if [[ "$GFLAG" == "-g" ]]; then
  # gzip mode: output is .gz; gzip level is only relevant here
  read d_time d_peak_mb < <(
    run_tm "decompress" spring -d -i "$SPRING_OUT" -o "$DEC_OUT" -g \
      --gzip-level "$OUTPUT_GZIP_LEVEL" -t "$THREADS" -w "$TMPDIR"
  )
else
  read d_time d_peak_mb < <(
    run_tm "decompress" spring -d -i "$SPRING_OUT" -o "$DEC_OUT" \
      -t "$THREADS" -w "$TMPDIR"
  )
fi
[[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }

# Decompressed size
dec_bytes="$(stat -c%s "$DEC_OUT")"
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"

# Verify integrity
stream() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -cd -- "$f"
  else
    cat -- "$f"
  fi
}

input_matches_decompressed=0
if cmp -s <(stream "$IN") <(stream "$DEC_OUT"); then
  input_matches_decompressed=1
fi
echo "input_matches_decompressed=$input_matches_decompressed"
# Metrics TSV
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