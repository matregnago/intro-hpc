#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=02:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=parsec


set -euo pipefail

echo "Date    : $(date)"
echo "Node(s) : $SLURM_JOB_NODELIST"
echo "JobID   : $SLURM_JOB_ID"

WORK_DIR="$SCRATCH/intro_hpc"
RESULTS_DIR="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$FINAL_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

trap 'cp -r "$WORK_DIR/$RESULTS_DIR/." "$FINAL_DIR/" 2>/dev/null || true' EXIT

export THREADS=16 GPUS=1
export DESIGN_FILE="scripts/doe/doe_parsec_dag.csv"

nixw nix develop .#parsec --impure --command env \
    RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
    THREADS="$THREADS" GPUS="$GPUS" \
    bash scripts/capture_parsec_dag.sh

echo "==> Pronto. Resultados em data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/runs/"
