#!/bin/bash -l
#SBATCH --job-name=n_size_tupi
#SBATCH --partition=tupi
#SBATCH --exclude=tupi1,tupi2
#SBATCH --time=11:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euxo pipefail

WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/local_results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_RESULTS_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
rsync -a \
  --exclude='/.git/' \
  --exclude='/.direnv/' \
  --exclude='/data/' \
  --exclude='/papers/' \
  --exclude='/slides/' \
  --exclude='/plots/' \
  --exclude='/context/' \
  --exclude='/result' \
  --exclude='/result-*' \
  --exclude='/*.out' \
  --exclude='/*.err' \
  "$SLURM_SUBMIT_DIR/" "$WORK_DIR/"
cd "$WORK_DIR"

copyback() {
    echo "$LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    if [ -d "$LOCAL_RESULTS_DIR" ]; then
        cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
    fi
}
trap copyback EXIT

export RESULTS_DIR="$LOCAL_RESULTS_DIR"
mkdir -p "$RESULTS_DIR"
cp doe_n_size_tupi.csv "$RESULTS_DIR/design.csv"
export DESIGN_FILE="$RESULTS_DIR/design.csv"

export THREADS=23
export GPUS=1

nixw nix develop .#starpu --command bash -c '
set -euo pipefail
bash scripts/run.sh starpu
'

nixw nix develop .#parsec --command bash -c '
set -euo pipefail
bash scripts/run.sh parsec
'
