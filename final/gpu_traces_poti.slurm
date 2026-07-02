#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=04:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=gpu_traces_poti

# Captura de rastros (StarVZ) na poti (i7-14700KF 20c, RTX 4070 12GB) para os
# dois kernels FP64 do estudo (dpotrf, dgetrf_nopiv) com b=500 (pico de E1),
# N=40000, THREADS=19 (reserva 1 core p/ dirigir a GPU), 1 GPU.
#
# DoE: final/doe_gpu_traces_poti_500.csv (8 runs: 2 algoritmos x 4
# runtime/scheduler x 1 rep). Roda nos dois runtimes com TRACE=1 + TRACE_FULL=1
# + TRACE_STATS=1:
#   - StarPU: FxT + BUS_STATS + WORKER_STATS (FxT ja e completo sob TRACE=1).
#   - PaRSEC: profile (.prof) + task_profiler PINS (states EXEC/RELEASE_DEPS/
#     ACTIVATE_CB no .prof) + device_show_statistics.
# TRACE_DAG=0 -- o .dot fica para um job dedicado (slurm/parsec_dag.slurm),
# pois o grapher polui o rastro cronometrado.
#
# Os modelos de desempenho do StarPU (dmda/dmdas) sao CALIBRADOS antes das
# medidas (final/gpu_doe_sweep.sh::calibrate_starpu, CALIB_PASSES=2, modelo
# isolado em $RESULTS_DIR/.starpu, congelado p/ as runs cronometradas com
# STARPU_CALIBRATE=0).
#
#   sbatch final/gpu_traces_poti.slurm
#
# Saida: data/gpu_traces_poti_<jobid>/{results.csv,design.csv,runs/.../,
#       .starpu/}  (+  <jobid>.{out,err} no $SLURM_SUBMIT_DIR)

set -euo pipefail

NODE=poti
THREADS=19
GPUS=1
DESIGN_SRC="final/doe_gpu_traces_poti_500.csv"

echo "==============================================================="
echo "Date        : $(date)"
echo "JobID       : $SLURM_JOB_ID   Partition: $SLURM_JOB_PARTITION"
echo "Node(s)     : $SLURM_JOB_NODELIST   THREADS=$THREADS GPUS=$GPUS"
echo "==============================================================="

WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_RESULTS_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

[[ -f "$DESIGN_SRC" ]] || { echo "DoE nao encontrado: $DESIGN_SRC" >&2; exit 1; }

copyback() {
    echo "$LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    [ -d "$LOCAL_RESULTS_DIR" ] && cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
}
trap copyback EXIT

export RESULTS_DIR="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
cp "$DESIGN_SRC" "${RESULTS_DIR}/design.csv"
export DESIGN_FILE="${RESULTS_DIR}/design.csv"

# Rastro completo + contadores nativos de movimentacao de dados. Sem DAG.
export TRACE=1 TRACE_FULL=1 TRACE_DAG=0 TRACE_STATS=1
export CALIB_PASSES=2
export THREADS GPUS

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> fase: ${runtime}  (nixw nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" \
        TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" TRACE_DAG="$TRACE_DAG" \
        TRACE_STATS="$TRACE_STATS" CALIB_PASSES="$CALIB_PASSES" \
        bash final/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! fase ${runtime} (.#${shell}) falhou; continuando" >&2
    fi
}

# StarPU primeiro (calibra dmda/dmdas), depois PaRSEC.
run_phase starpu starpu
run_phase parsec parsec

echo "==> rastros (poti) concluidos."
