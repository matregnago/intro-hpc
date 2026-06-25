#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=02:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=traces

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

cp "scripts/doe/doe_full_trace.csv" "${RESULTS_DIR}/design.csv"
export DESIGN_FILE="${RESULTS_DIR}/design.csv"

# Captura completa do trace. TRACE_DAG=0: a DAG sai do job pequeno dedicado
# (slurm/parsec_dag.slurm); em n=19200 ela seria irrenderavel e o grapher
# poluiria o tempo do trace. Mantemos TRACE_FULL=1 pelo PINS task_profiler,
# que enriquece o .prof do PaRSEC para o StarVZ.
export TRACE=1
export TRACE_FULL=1
export TRACE_DAG=0
export THREADS=16
export GPUS=1

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> phase: ${runtime}  (nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" \
        TRACE_DAG="$TRACE_DAG" \
        bash scripts/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! phase ${runtime} (.#${shell}) failed; continuing" >&2
    fi
}

run_phase starpu starpu
run_phase parsec parsec
