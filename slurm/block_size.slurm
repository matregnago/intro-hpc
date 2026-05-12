#!/bin/bash -l
#SBATCH --partition=cei
#SBATCH --time=23:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=chameleon-runtimes

set -euo pipefail

echo "==============================================================="
echo "Date        : $(date)"
echo "SubmitHost  : $(hostname)"
echo "SubmitDir   : $SLURM_SUBMIT_DIR"
echo "JobID       : $SLURM_JOB_ID"
echo "Node(s)     : $SLURM_JOB_NODELIST"
echo "==============================================================="


WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_RESULTS_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

copyback() {
    echo "$LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    if [ -d "$LOCAL_RESULTS_DIR" ]; then
        cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
    fi
}
trap copyback EXIT

export RESULTS_DIR="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"

DESIGN_FILE="${SLURM_SUBMIT_DIR}/doe_block_size.csv"

cp "$DESIGN_FILE" "${RESULTS_DIR}/doe_block_size.csv"
export DESIGN_FILE="${RESULTS_DIR}/doe_block_size.csv"


# StarPU
nixw nix develop .#starpu --command bash -c '
set -euo pipefail
bash scripts/run.sh starpu
'

# PaRSEC
nixw nix develop .#parsec --command bash -c '
set -euo pipefail
bash scripts/run.sh parsec
'
