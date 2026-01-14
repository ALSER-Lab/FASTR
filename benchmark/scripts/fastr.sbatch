#!/usr/bin/env bash
#SBATCH --account=csc1685s267
#SBATCH --job-name=fastr_r1
#SBATCH --output=/home/users/aadeniyi6/logs/fastr_r1_%A_%a.out
#SBATCH --error=/home/users/aadeniyi6/logs/fastr_r1_%A_%a.err
#SBATCH --partition=qDEV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aadeniyi6@student.gsu.edu

set -euo pipefail

mkdir -p /home/users/aadeniyi6/logs

# Inputs
R1="/home/users/aadeniyi6/sra_downloads/fastq/SRR36648629/SRR36648629_1.fastq"
OUTROOT="/home/users/aadeniyi6/bench_results/fastr/SRR36648629_R1"

# FASTR scripts
TO_FASTR="/home/users/aadeniyi6/FASTR/src/to_fastr.py"
TO_FASTQ="/home/users/aadeniyi6/FASTR/src/to_fastq.py"

# Wrapper
FASTR_WRAPPER="/home/users/aadeniyi6/scripts/fastr.sh"

# Params
MODE="${SLURM_ARRAY_TASK_ID}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
TMPDIR="/tmp"

QUAL_SCALE="log"
SEQ_TYPE="illumina_sra"
PHRED_ALPHA="phred94"

# Conda
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/users/aadeniyi6/fastq_bench

# Sanity checks
[[ -x "$FASTR_WRAPPER" ]] || { echo "ERROR: wrapper not executable: $FASTR_WRAPPER"; exit 1; }
[[ -f "$TO_FASTR" ]] || { echo "ERROR: missing: $TO_FASTR"; exit 1; }
[[ -f "$TO_FASTQ" ]] || { echo "ERROR: missing: $TO_FASTQ"; exit 1; }
[[ -d "$TMPDIR" && -w "$TMPDIR" ]] || { echo "ERROR: TMPDIR must exist+writable: $TMPDIR"; exit 1; }

# Confirm deps early (to_fastq needs numba)
python -c "import numpy, numba; print('deps ok:', numpy.__version__)"  # will fail-fast if missing

mkdir -p "$(dirname "$OUTROOT")"

echo "Host=$(hostname) Start=$(date -Is) MODE=$MODE THREADS=$THREADS"
echo "R1=$R1"
echo "OUTROOT=$OUTROOT"
echo "TO_FASTR=$TO_FASTR"
echo "TO_FASTQ=$TO_FASTQ"
echo "QUAL_SCALE=$QUAL_SCALE SEQ_TYPE=$SEQ_TYPE PHRED_ALPHA=$PHRED_ALPHA"
echo

bash "$FASTR_WRAPPER" \
  --input "$R1" \
  --mode "$MODE" \
  --threads "$THREADS" \
  --tmpdir "$TMPDIR" \
  --to-fastr "$TO_FASTR" \
  --to-fastq "$TO_FASTQ" \
  --qual_scale "$QUAL_SCALE" \
  --seq_type "$SEQ_TYPE" \
  --phred_alpha "$PHRED_ALPHA" \
  --decompress-to-null \
  --no-verify \
  --out-prefix "${OUTROOT}_mode${MODE}"

echo "Done $(date -Is)"
