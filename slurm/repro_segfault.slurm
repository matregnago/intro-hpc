#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --gres=gpu:1
#SBATCH --time=00:15:00
#SBATCH --job-name=reprosegf
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.out

set -uo pipefail
echo "Date    : $(date)"
echo "Node(s) : $SLURM_JOB_NODELIST"
echo "JobID   : $SLURM_JOB_ID"

cd "$SLURM_SUBMIT_DIR"

nixw nix develop .#parsec --impure --command bash scripts/repro_segfault_inner.sh
echo "==> repro done"
