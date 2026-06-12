#!/bin/bash -l
#SBATCH --partition=cidia
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

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

# StarPU
nixw nix develop .#starpu-cuda --impure --command bash -c '
set -euo pipefail
date
'
