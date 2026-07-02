#!/bin/bash -l
#SBATCH --partition=tupi
#SBATCH --time=04:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --job-name=gpu_trace_ideal_tupi
#SBATCH --nodelist=tupi3,tupi4,tupi5,tupi6

# E4-ideal -- COLETA COMPLETA de rastro StarVZ na GPU, particao TUPI
# (i9-14900KF 24c, RTX 4090 24GB). THREADS=23 (reserva 1 core p/ dirigir a GPU).
#
# Roda final/doe_gpu_trace_ideal_tupi.csv (16 execucoes: 4 sched x 2 algos x 2
# precisoes, 1 rep) nos dois runtimes, no TILE IDEAL de cada (precisao,algoritmo)
# -- MESMO design/tiles do poti (os tiles vem do E1b poti; a tupi ainda nao tem
# run limpo). Captura TUDO p/ o StarVZ:
#   TRACE=1       -> StarPU FxT (.fxt) / PaRSEC .prof
#   TRACE_FULL=1  -> PINS task_profiler no PaRSEC (enriquece o .prof) p/ os
#                    converters dag/tasks/states.parquet
#   TRACE_DAG=1   -> DAG .dot do PaRSEC (CHAMELEON_PARSEC_DOT) p/ dag.parquet +
#                    reconstrucao do Ready. ATENCAO: nesses N a DAG e enorme
#                    (so usavel via pipeline parquet, nao renderavel como grafo)
#                    e o grapher adiciona overhead ao tempo medido.
#   TRACE_STATS=1 -> contadores nativos de data-movement no log.
# O StarPU roda com 4 workers CUDA (STARPU_NWORKER_PER_CUDA=4, default do
# gpu_doe_sweep.sh) p/ casar com os 4 exec streams do PaRSEC -- assim os paineis
# Espaco-Tempo do StarVZ tem o mesmo numero de lanes de GPU nos dois runtimes.
# Os modelos do StarPU (dmda/dmdas) sao CALIBRADOS antes das medidas
# (final/gpu_doe_sweep.sh::calibrate_starpu, modelo isolado em
# $RESULTS_DIR/.starpu, congelado p/ as runs cronometradas).
#
# OBS: o run anterior da tupi (797982) teve a GPU aparentemente disputada
# (StarPU OOM na calibracao de banda + GFlops do PaRSEC abaixo do 4070); confira
# que a GPU do no esta livre antes de submeter.
#
#   sbatch final/gpu_trace_ideal_tupi.slurm

set -euo pipefail

NODE=tupi
THREADS=23
GPUS=1
DESIGN_SRC="final/doe_gpu_trace_ideal_${NODE}.csv"

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

# Coleta COMPLETA de rastro (trace + PINS + DAG + stats nativos).
export TRACE=1 TRACE_FULL=1 TRACE_DAG=1 TRACE_STATS=1
export THREADS GPUS

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> fase: ${runtime}  (nixw nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" \
        TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" TRACE_DAG="$TRACE_DAG" \
        TRACE_STATS="$TRACE_STATS" \
        bash final/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! fase ${runtime} (.#${shell}) falhou; continuando" >&2
    fi
}

# StarPU primeiro (calibra dmda/dmdas), depois PaRSEC.
run_phase starpu starpu
run_phase parsec parsec

echo "==> E4-ideal (tupi) concluido."
