#!/bin/bash

# Roda NO NO de computacao (tupi), dentro da sua alocacao interativa
# (`sbatch -p tupi --gres=gpu:1 --cpus-per-task=24 --wrap 'sleep infinity'`).
# Equivalente interativo do slurm/gpu_tile_tupi.slurm mas para o EXPERIMENTO
# E2b (escalonamento vs tamanho do problema N, no pico de tile de E1 da tupi).
#
# DoE: final/doe_gpu_n_size_tupi.csv (240 runs: 2 alg FP64 x 4 runtime/scheduler
# x 10 Ns x 3 reps, todos com b=1000 -- o pico de tile FP64 achado em E1 / job
# gpu_tile_tupi_20260701_095748: StarPU pica em b=1000). Sem traces (so
# time/gflops) -- TRACE=0. Modelos dmda/dmdas do StarPU sao CALIBRADOS antes das
# medidas (CALIB_PASSES=2, modelo isolado em $RESULTS_DIR/.starpu, congelado p/
# as runs cronometradas com STARPU_CALIBRATE=0).
#
# MODO EXTEND: para complementar um job anterior (Ns 5000-40000 ja medidos em
# data/gpu_n_size_tupi_20260701_115305), gere um CSV so com os Ns novos e aponte
# DOE_CSV para ele. O calibrate_starpu so calibra (precision,algorithm,n,b)
# ainda nao vistos, entao roda so o necessario.
#   OUT=final/doe_gpu_n_size_tupi_extend.csv N_SIZES="50000 60000" \
#     Rscript final/doe_gpu_n_size_tupi.r
#   DOE_CSV=final/doe_gpu_n_size_tupi_extend.csv \
#     JOBTAG=extend_$(date +%Y%m%d_%H%M%S) \
#     bash scripts/pcad/collect_n_size_tupi_node.sh
# Depois concatene results.csv do job extend ao do job original (mesmo schema).
#
# Normalmente disparado pelo front-end (scripts/pcad/collect_n_size_tupi_pcad.sh),
# mas roda direto se voce ja esta no no:
#   bash scripts/pcad/collect_n_size_tupi_node.sh
#
# Env: SUBMIT_DIR (repo no home, default ~/intro_hpc), SCRATCH, JOBTAG,
#      THREADS (override, default 23), GPUS (1).

set -euo pipefail

SUBMIT_DIR="${SUBMIT_DIR:-$HOME/intro_hpc}"
JOBTAG="${JOBTAG:-$(date +%Y%m%d_%H%M%S)}"

HOST="$(hostname)"
case "$HOST" in
    *tupi*) NODE_TAG="tupi"; DEFAULT_THREADS=23 ;;
    *)    echo "!! no desconhecido: $HOST (esperado tupi)" >&2; exit 1 ;;
esac
DOE_CSV="${DOE_CSV:-final/doe_gpu_n_size_tupi.csv}"

if [[ -z "${SCRATCH:-}" ]]; then
    SCRATCH="$HOME/scratch"
    echo "!! \$SCRATCH nao definido; usando $SCRATCH (convencao do nixw)." >&2
    echo "   Prefira disparar via srun --jobid (herda o \$SCRATCH da alocacao)." >&2
fi

WORK_DIR="$SCRATCH/intro_hpc"
RESULTS_REL="data/gpu_n_size_${NODE_TAG}_${JOBTAG}"
LOCAL_RESULTS_DIR="$WORK_DIR/$RESULTS_REL"
FINAL_RESULTS_DIR="$SUBMIT_DIR/$RESULTS_REL"

echo "==============================================================="
echo "Node       : $HOST  (tag=$NODE_TAG)"
echo "SubmitDir  : $SUBMIT_DIR"
echo "WorkDir    : $WORK_DIR   (\$SCRATCH=$SCRATCH)"
echo "JobTag     : $JOBTAG"
echo "DoE        : $DOE_CSV"
echo "==============================================================="

[[ -f "$SUBMIT_DIR/flake.nix" ]] || { echo "repo nao encontrado em $SUBMIT_DIR" >&2; exit 1; }

mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
# Espelha o repo no scratch. Exclui os dados pesados e o .git/.direnv.
rsync -a --delete \
    --exclude='.git/' --exclude='.direnv/' --exclude='.claude/' \
    --exclude='data/' --exclude='plots/' --exclude='context/' \
    "$SUBMIT_DIR/" "$WORK_DIR/"
cd "$WORK_DIR"

copyback() {
    echo ">> copyback: $LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    if [[ -d "$LOCAL_RESULTS_DIR" ]]; then
        mkdir -p "$FINAL_RESULTS_DIR"
        cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
    fi
}
trap copyback EXIT

export RESULTS_DIR="$RESULTS_REL"   # relativo ao WORK_DIR
export THREADS="${THREADS:-$DEFAULT_THREADS}"
export GPUS="${GPUS:-1}"

# DoE de N-size para a tupi (FP64, b=1000, 10 Ns x 2 alg x 4 sched x 3 reps = 240 runs).
cp "$DOE_CSV" "${RESULTS_REL}/design.csv"
export DESIGN_FILE="${RESULTS_REL}/design.csv"

# Sem traces nesta fase -- so medimos time/gflops.
export TRACE=0 TRACE_FULL=0 TRACE_DAG=0 TRACE_STATS=0
export CALIB_PASSES=2

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> phase: ${runtime}  (nixw nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" \
        TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" TRACE_DAG="$TRACE_DAG" \
        TRACE_STATS="$TRACE_STATS" CALIB_PASSES="$CALIB_PASSES" \
        bash final/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! phase ${runtime} (.#${shell}) failed; continuing" >&2
    fi
}

# StarPU primeiro (calibra dmda/dmdas), depois PaRSEC.
run_phase starpu starpu
run_phase parsec parsec

echo "==> Resultados em $FINAL_RESULTS_DIR  (no \$HOME, compartilhado com o front-end)"
