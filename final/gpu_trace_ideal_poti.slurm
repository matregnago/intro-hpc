#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=04:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=gpu_trace_ideal_poti

# E4-ideal -- COLETA COMPLETA de rastro StarVZ na GPU, particao POTI
# (i7-14700KF 20c, RTX 4070 12GB). THREADS=19 (reserva 1 core p/ dirigir a GPU).
#
# Roda final/doe_gpu_trace_ideal_poti.csv (16 execucoes: 4 sched x 2 algos x 2
# precisoes, 1 rep) nos dois runtimes, no TILE IDEAL de cada (precisao,algoritmo)
# medido no E1b. Captura TUDO p/ o StarVZ:
#   TRACE=1       -> StarPU FxT (.fxt) / PaRSEC .prof
#   TRACE_FULL=1  -> PINS task_profiler no PaRSEC (enriquece o .prof: EXEC/
#                    RELEASE_DEPS/etc) p/ os converters dag/tasks/states.parquet
#   TRACE_DAG=1   -> DAG .dot do PaRSEC (CHAMELEON_PARSEC_DOT) p/ dag.parquet +
#                    reconstrucao do Ready. ATENCAO: nesses N a DAG e enorme
#                    (so usavel via pipeline parquet, nao renderavel como grafo)
#                    e o grapher adiciona overhead ao tempo medido -- aqui o foco
#                    e a estrutura/StarVZ, nao o tempo de parede.
#   TRACE_STATS=1 -> contadores nativos de data-movement no log (StarPU BUS/
#                    WORKER_STATS, PaRSEC device_show_statistics).
# O StarPU roda com 4 workers CUDA (STARPU_NWORKER_PER_CUDA=4, default do
# gpu_doe_sweep.sh) p/ casar com os 4 exec streams do PaRSEC -- assim os paineis
# Espaco-Tempo do StarVZ tem o mesmo numero de lanes de GPU nos dois runtimes.
# Os modelos do StarPU (dmda/dmdas) sao CALIBRADOS por (precisao,algoritmo,n,b)
# antes das medidas (final/gpu_doe_sweep.sh::calibrate_starpu, modelo isolado em
# $RESULTS_DIR/.starpu, congelado p/ as runs cronometradas).
#
#   sbatch final/gpu_trace_ideal_poti.slurm

set -euo pipefail

NODE=poti
THREADS=19
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

echo "==> E4-ideal (poti) concluido."
