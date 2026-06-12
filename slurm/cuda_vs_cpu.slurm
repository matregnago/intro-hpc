#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=08:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=chameleon-cuda-vscpu
# poti já expõe o RTX 4070 ao job; se algum dia a partição exigir alocação
# explícita de GPU, descomente a linha abaixo:
##SBATCH --gres=gpu:1

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

# Sweep shared by all four phases (overridable here).
# FP32 (`s`): the RTX 4070 has ~29 TFlops FP32 vs ~1.7 TFlops on the 16-thread
# CPU, so the GPU's acceleration is clearly visible here (sgemm ~4.4x vs CPU-16
# in docs/parsec_cuda_results.md). In FP64 (`d`) the card is hardware-capped
# below the CPU (~1x), which hides whether the GPU is doing anything — so this
# "is the GPU accelerating?" test runs in `s`. sgemm is the cleanest signal
# (all trailing-update); spotrf/sgeqrf accelerate less because their panel +
# vertical-merge critical path stays on CPU by design (matches StarPU).
# Switch to "dpotrf dgeqrf" when you want the double-precision numbers to report.
export OPS="sgemm spotrf sgeqrf"
# n estendido p/ ver escalabilidade; b=2000 favorece os tiles do cuBLAS na GPU.
# Todos os n sao multiplos de 2000 -> tiling limpo. CPU em n=48000 ~168s/rep de
# sgemm (O(n^3)); GPU continua barata. Wall estimado ~1,5-2h (limite 08:00:00).
export SWEEP_PAIRS="8000:2000 16000:2000 24000:2000 32000:2000 40000:2000 48000:2000"
export REPS=3
export CPU_THREADS=16
export GPU_THREADS=8
# VERBOSE=1 confirms GPU tile dispatch in the logs (noisy at large n).
export VERBOSE=0

# Run one (runtime, mode) phase inside its dev shell. Phases are independent:
# a failure in one (e.g. a missing build) is logged but does not abort the
# rest, so a single submission still produces whatever results did succeed.
#   --impure is required: chameleon.nix reads the local source tree via
#   builtins.getEnv "PWD", so the eval must see the rsync'd WORK_DIR.
run_phase() {
    local shell="$1" runtime="$2" mode="$3"
    echo ">>> phase: ${runtime} ${mode}  (nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command bash -c "
        set -euo pipefail
        bash scripts/run_cuda_vscpu.sh ${runtime} ${mode}
    "; then
        echo "!!! phase ${runtime} ${mode} (.#${shell}) failed; continuing" >&2
    fi
}

# StarPU
run_phase starpu      starpu cpu   # CPU-only StarPU
run_phase starpu-cuda starpu gpu   # StarPU + CUDA

# PaRSEC
run_phase parsec-cpu  parsec cpu   # CPU-only PaRSEC
run_phase parsec      parsec gpu   # PaRSEC + CUDA (cuBLAS/cuSOLVER)
