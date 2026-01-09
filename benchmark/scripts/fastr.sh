#!/usr/bin/env bash
set -euo pipefail

# FASTR benchmark wrapper (single FASTQ)
# Compress:
#   python to_fastr.py IN.fastq OUT.fastr --mode 0..3 --qual_scale x --seq_type illumina_sra --workers N --phred_alpha phred94
# Decompress:
#   python to_fastq.py OUT.fastr OUT.fastq --mode M --num_workers N --phred_alpha phred94 [--headers_file ...]  (only mode=3)

IN=""
THREADS="16"
TMPDIR="/tmp"
OUT_PREFIX="fastr-output"
MODE=""

TO_FASTR=""   # /path/to/to_fastr.py
TO_FASTQ=""   # /path/to/to_fastq.py

QUAL_SCALE="x"
SEQ_TYPE="illumina_sra"
PHRED_ALPHA="phred94"

DECOMPRESS_TO_NULL=0
VERIFY=0

TOOL="fastr"

THREAD_FLAG_COMPRESS="--workers"
THREAD_FLAG_DECOMPRESS="--num_workers"
MT_COMPRESS=1
MT_DECOMPRESS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) IN="$2"; shift 2;;
    --threads|--workers) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;

    --to-fastr|--fastr-main) TO_FASTR="$2"; shift 2;;   # keep --fastr-main as alias
    --to-fastq) TO_FASTQ="$2"; shift 2;;

    --qual_scale|--qual-scale) QUAL_SCALE="$2"; shift 2;;
    --seq_type|--seq-type) SEQ_TYPE="$2"; shift 2;;
    --phred_alpha|--phred-alpha) PHRED_ALPHA="$2"; shift 2;;

    --decompress-to-null) DECOMPRESS_TO_NULL=1; shift 1;;
    --verify) VERIFY=1; shift 1;;
    --no-verify) VERIFY=0; shift 1;;

    -h|--help)
      cat <<EOF
Usage:
  $0 --input <file.fastq[.gz]> --mode 0|1|2|3 --to-fastr /path/to/to_fastr.py [--to-fastq /path/to/to_fastq.py]
     [--threads N] [--tmpdir DIR] [--out-prefix PREFIX]
     [--qual_scale x] [--seq_type illumina_sra] [--phred_alpha phred94]
     [--decompress-to-null] [--verify|--no-verify]

Notes:
  - If --to-fastq is not provided, we assume it lives next to to_fastr.py as: ../to_fastq.py
  - Only mode=3 requires --headers_file on decode; this wrapper auto-detects the headers file.
EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1;;
  esac
done

[[ -n "$IN" ]] || { echo "ERROR: --input is required"; exit 1; }
[[ -n "$MODE" ]] || { echo "ERROR: --mode is required"; exit 1; }
[[ "$MODE" =~ ^[0-3]$ ]] || { echo "ERROR: --mode must be 0,1,2,3"; exit 1; }

[[ -n "$TO_FASTR" ]] || { echo "ERROR: --to-fastr is required"; exit 1; }
[[ -f "$TO_FASTR" ]] || { echo "ERROR: to_fastr.py not found: $TO_FASTR"; exit 1; }

if [[ -z "$TO_FASTQ" ]]; then
  TO_FASTQ="$(dirname "$TO_FASTR")/to_fastq.py"
fi
[[ -f "$TO_FASTQ" ]] || { echo "ERROR: to_fastq.py not found: $TO_FASTQ"; exit 1; }

OUT_PREFIX="${OUT_PREFIX%/}"
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads/--workers must be an integer"; exit 1; }
[[ -d "$TMPDIR" && -w "$TMPDIR" ]] || { echo "ERROR: --tmpdir must exist and be writable: $TMPDIR"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "ERROR: python not found in PATH"; exit 1; }

OUT_DIR="$(dirname "$OUT_PREFIX")"
BASE_NAME="$(basename "$OUT_PREFIX")"
[[ "$OUT_DIR" != "." ]] && mkdir -p "$OUT_DIR"
LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGFILE="$LOGDIR/${TOOL}_${BASE_NAME}_${RUN_ID}.log"
METRICS_TSV="$LOGDIR/metrics.tsv"

exec > >(tee -a "$LOGFILE") 2>&1

FASTR_OUT="${OUT_PREFIX}.mode${MODE}.fastr"
DEC_OUT="${OUT_PREFIX}.dec.mode${MODE}.fastq"

echo "=== FASTR benchmark ==="
echo "Start: $(date -Is)"
echo "tool=$TOOL"
echo "run_id=$RUN_ID"
echo "input=$IN"
echo "mode=$MODE"
echo "workers=$THREADS"
echo "qual_scale=$QUAL_SCALE"
echo "seq_type=$SEQ_TYPE"
echo "phred_alpha=$PHRED_ALPHA"
echo "to_fastr=$TO_FASTR"
echo "to_fastq=$TO_FASTQ"
echo "tmpdir=$TMPDIR"
echo "out_prefix=$OUT_PREFIX"
echo "decompress_to_null=$DECOMPRESS_TO_NULL"
echo "verify=$VERIFY"
echo "logfile=$LOGFILE"
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

stream() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then gzip -cd -- "$f"; else cat -- "$f"; fi
}

echo "Computing raw input bytes (uncompressed)..."
in_bytes="$(stream "$IN" | wc -c)"
echo "input_raw_bytes=$in_bytes"
echo


# Compress: FASTQ -> FASTR

rm -f "$FASTR_OUT"
read c_time c_peak_mb < <(
  run_tm "compress" python "$TO_FASTR" "$IN" "$FASTR_OUT" \
    --mode "$MODE" \
    --qual_scale "$QUAL_SCALE" \
    --seq_type "$SEQ_TYPE" \
    --workers "$THREADS" \
    --phred_alpha "$PHRED_ALPHA"
)
[[ -s "$FASTR_OUT" ]] || { echo "ERROR: compression did not produce $FASTR_OUT"; exit 1; }

out_bytes="$(stat -c%s "$FASTR_OUT")"
comp_size_mb="$(awk -v b="$out_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
comp_ratio="$(awk -v inb="$in_bytes" -v outb="$out_bytes" 'BEGIN{printf "%.6f", inb/outb}')"


# Mode 3 header file detection (decode requires it)

HEADERS_FILE=""
if [[ "$MODE" == "3" ]]; then
  # common candidates seen in your example
  cand1="${OUT_PREFIX}_headers.txt"
  cand2="${FASTR_OUT%.fastr}_headers.txt"
  cand3="${OUT_PREFIX}.headers.txt"

  if [[ -s "$cand1" ]]; then
    HEADERS_FILE="$cand1"
  elif [[ -s "$cand2" ]]; then
    HEADERS_FILE="$cand2"
  elif [[ -s "$cand3" ]]; then
    HEADERS_FILE="$cand3"
  else
    echo "ERROR: mode=3 requires a headers file for decoding, but none found."
    echo "Tried:"
    echo "  $cand1"
    echo "  $cand2"
    echo "  $cand3"
    echo "List output dir:"
    ls -lah "$(dirname "$OUT_PREFIX")" || true
    exit 1
  fi
  echo "headers_file=$HEADERS_FILE"
fi


# Decompress: FASTR -> FASTQ (or FIFO)

input_matches_decompressed=0

decode_cmd_base=(python "$TO_FASTQ" "$FASTR_OUT")
decode_flags=(--mode "$MODE" --num_workers "$THREADS" --phred_alpha "$PHRED_ALPHA")
if [[ "$MODE" == "3" ]]; then
  decode_flags+=(--headers_file "$HEADERS_FILE")
fi

if [[ $DECOMPRESS_TO_NULL -eq 1 ]]; then
  FIFO="$TMPDIR/fastr_fifo_${RUN_ID}_$$"
  BYTES_FILE="$TMPDIR/fastr_dec_bytes_${RUN_ID}_$$.txt"
  rm -f "$FIFO" "$BYTES_FILE"
  mkfifo "$FIFO"

  (wc -c < "$FIFO" > "$BYTES_FILE") &
  WC_PID=$!

  read d_time d_peak_mb < <(
    run_tm "decompress" "${decode_cmd_base[@]}" "$FIFO" "${decode_flags[@]}"
  )

  wait "$WC_PID" || true
  dec_bytes="$(tr -d '[:space:]' < "$BYTES_FILE")"
  rm -f "$BYTES_FILE" "$FIFO"
else
  rm -f "$DEC_OUT"
  read d_time d_peak_mb < <(
    run_tm "decompress" "${decode_cmd_base[@]}" "$DEC_OUT" "${decode_flags[@]}"
  )
  [[ -s "$DEC_OUT" ]] || { echo "ERROR: decompression did not produce $DEC_OUT"; exit 1; }
  dec_bytes="$(stat -c%s "$DEC_OUT")"

  if [[ $VERIFY -eq 1 ]]; then
    if cmp -s <(stream "$IN") "$DEC_OUT"; then input_matches_decompressed=1; fi
  fi
fi

dec_size_mb="$(awk -v b="$dec_bytes" 'BEGIN{printf "%.3f", b/1048576.0}')"
echo "input_matches_decompressed=$input_matches_decompressed"


# Metrics TSV

if [[ ! -f "$METRICS_TSV" ]]; then
  printf "tool\trun_id\tthreads_requested\tthread_flag_compress\tthread_flag_decompress\tmt_compress\tmt_decompress\tcompression_time_s\tcompressed_size_mb\tcompression_ratio\tpeak_mem_compress_mb\tdecompression_time_s\tdecompressed_size_mb\tpeak_mem_decompress_mb\tinput_matches_decompressed\n" \
    >> "$METRICS_TSV"
fi

printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.6f\t%s\t%.6f\t%.3f\t%.6f\t%s\t%.3f\t%d\n" \
  "${TOOL}_mode${MODE}" "$RUN_ID" "$THREADS" "$THREAD_FLAG_COMPRESS" "$THREAD_FLAG_DECOMPRESS" \
  "$MT_COMPRESS" "$MT_DECOMPRESS" \
  "$c_time" "$comp_size_mb" "$comp_ratio" "$c_peak_mb" \
  "$d_time" "$dec_size_mb" "$d_peak_mb" \
  "$input_matches_decompressed" \
  | tee -a "$METRICS_TSV"

echo
echo "Done: $(date -Is)"
echo "Log saved to: $LOGFILE"
echo "Metrics appended to: $METRICS_TSV"
