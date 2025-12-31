#!/usr/bin/env bash
set -euo pipefail

# Wrapper for Info-ZIP benchmarking.
# Runs compress -> decompress cycle on a single file.
#
# Note: Maps --fast/--best to zip -1/-9 
#
# Usage:
#   ./zip.sh --input file.fastq --fast --out-prefix /path/to/results/run1

#
# Configs and Defaults
#
IN=""
MODE="fast"           # fast (-1) or best (-9)
TMPDIR="/tmp"
OUT_PREFIX="zip-output"
TOOL="zip"

# zip/unzip is single-threaded
THREADS="1"
MT_COMPRESS=0
MT_DECOMPRESS=0
THREAD_FLAG_COMPRESS="NA"
THREAD_FLAG_DECOMPRESS="NA"

#
# Parse named args
#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --fast) MODE="fast"; shift 1;;
    --best) MODE="best"; shift 1;;
    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file> [--fast|--best] [--tmpdir DIR] [--out-prefix PREFIX]

Notes:
  - Default is --fast (zip -1).
  - Use OUT_PREFIX to encode the run identifier (e.g., .../run1_fast, .../run2_best).
EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1;;
  esac
done

[[ -n "$IN" ]] || { echo "ERROR: --input is required"; exit 1; }
OUT_PREFIX="${OUT_PREFIX%/}"

# Map mode to actual zip flags
ZIP_LEVEL="1"
if [[ "$MODE" == "best" ]]; then
  ZIP_LEVEL="9"
fi

# --- Sanity Checks ---
# We need standard zip/unzip, usually found in the conda env or system path.
command -v zip >/dev/null 2>&1 || { echo "ERROR: zip not found in PATH (activate conda env)."; exit 1; }
command -v unzip >/dev/null 2>&1 || { echo "ERROR: unzip not found in PATH (install unzip in the same conda env)."; exit 1; }
command -v zipinfo >/dev/null 2>&1 || { echo "ERROR: zipinfo not found in PATH (usually comes with unzip/zip)."; exit 1; }

# TMPDIR checks
[[ -d "$TMPDIR" ]] || { echo "ERROR: --tmpdir does not exist: $TMPDIR"; exit 1; }
[[ -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir not writable: $TMPDIR"; exit 1; }

# logging
OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"

LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"

# Pipe all stdout/stderr to logfile and console
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== zip benchmark (single file) ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "mode=$MODE (zip -$ZIP_LEVEL)"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "logfile=$LOGFILE"
echo

ZIP_OUT="${OUT_PREFIX}.zip"
DEC_OUT="${OUT_PREFIX}.dec"

echo "zip_out=$ZIP_OUT"
echo "dec_out=$DEC_OUT"
echo

# Helper to capture /usr/bin/time metrics without parsing clutter
run_tm() {
  local label="$1"; shift
  local tf
  tf="$(mktemp -p "$TMPDIR" "time_${label}_XXXX.txt")"

  { echo; echo "---- $label ----"; echo "CMD: $*"; echo "Start($label): $(date -Is)"; } >&2

  set +e
  # Using system time binary for consistent memory reporting (%M is max RSS in KB)
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

#
# Compression
#
rm -f "$ZIP_OUT"
# -j junk paths (store just filename), -q quiet
read c_time c_peak_mb < <(run_tm "compress" zip -q "-$ZIP_LEVEL" -j "$ZIP_OUT" "$IN")
[[ -s "$ZIP_OUT" ]] || { echo "ERROR: compression did not produce $ZIP_OUT"; exit 1; }

# calculate ratio
in_bytes="$(stat -c%s "$IN")"
out_bytes="$(stat -c%s "$ZIP_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"

#
# Decompression
#
rm -f "$DEC_OUT"
# Need the internal filename to extract specifically that file
member="$(zipinfo -1 "$ZIP_OUT" | head -n 1)"
[[ -n "$member" ]] || { echo "ERROR: could not determine zip member name"; exit 1; }
echo "zip_member=$member"

# Using 'unzip -p' > file to avoid interactive prompts or path issues
read d_time d_peak_mb < <(
  run_tm "decompress" bash -c 'unzip -p -- "$1" "$2" > "$3"' unzip_pipe "$ZIP_OUT" "$member" "$DEC_OUT"
)
[[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }

dec_bytes="$(stat -c%s "$DEC_OUT")"
dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"

# verify integrity
input_matches_decompressed=0
if cmp -s -- "$IN" "$DEC_OUT"; then
  input_matches_decompressed=1
fi
echo "input_matches_decompressed=$input_matches_decompressed"

#
# Metrics TSV save
#
if [[ ! -f "$METRICS_TSV" ]]; then
# Header
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
