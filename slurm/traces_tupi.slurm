#!/bin/bash -l
#SBATCH --job-name=traces_tupi
#SBATCH --partition=tupi
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# E4 -- captura de rastros StarVZ na tupi (RTX 4090 24GB) para os 2 kernels
# FP64 do estudo (dpotrf, dgetrf_nopiv) com b=1000 (pico de E1), N=40000,
# THREADS=23, 1 GPU. DoE: doe_traces_tupi.csv (8 runs: 2 algoritmos x 4
# runtime/scheduler x 1 rep).
#
# Mesmas flags de traces_poti.slurm: TRACE=1 TRACE_FULL=1 TRACE_DAG=0
# TRACE_STATS=1, CALIB_PASSES=2.
#
#   sbatch slurm/traces_tupi.slurm
# Saida: data/traces_tupi_<jobid>/{results.csv,design.csv,runs/.../,.starpu/}

set -euxo pipefail

WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/local_results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
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

export RESULTS_DIR="$LOCAL_RESULTS_DIR"
mkdir -p "$RESULTS_DIR"
cp doe_traces_tupi.csv "$RESULTS_DIR/design.csv"
export DESIGN_FILE="$RESULTS_DIR/design.csv"

export THREADS=23
export GPUS=1
export TRACE=1 TRACE_FULL=1 TRACE_DAG=0 TRACE_STATS=1
export CALIB_PASSES=2

# StarPU (calibra dmda/dmdas antes das medidas -- ver scripts/run.sh)
nixw nix develop .#starpu --impure --command bash -c '
set -euo pipefail
bash scripts/run.sh starpu
'

# PaRSEC
nixw nix develop .#parsec --impure --command bash -c '
set -euo pipefail
bash scripts/run.sh parsec
'
