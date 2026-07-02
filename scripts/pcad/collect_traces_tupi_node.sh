#!/bin/bash

# Roda NO NO de computacao (tupi), dentro da sua alocacao interativa
# (`sbatch -p tupi --gres=gpu:1 --cpus-per-task=24 --wrap 'sleep infinity'`).
# Espelho tupi do collect_traces_poti_node.sh: rsync do repo -> $SCRATCH, roda
# as duas fases (starpu/parsec) do gpu_doe_sweep.sh la (I/O rapido), e copia os
# resultados de volta pro $HOME via trap EXIT (copyback).
#
# DoE: final/doe_gpu_traces_tupi_1000.csv (8 runs: 2 algoritmos FP64 x 4
# runtime/scheduler x 1 rep, N=40000, b=1000). N e o maior do sweep de N-size
# da tupi (E2b, ainda subindo) e b=1000 e o pico de tile do StarPU em E1 --
# mesmo b usado no sweep de N-size da tupi (PaRSEC pica em b=2000, mas mantemos
# b uniforme para comparar timelines com a mesma granularidade de tile).
# Captura rastro completo (TRACE=1 + TRACE_FULL=1 + TRACE_STATS=1):
#   - StarPU: FxT + BUS_STATS + WORKER_STATS. FxT ja e completo sob TRACE=1.
#   - PaRSEC: profile (.prof) + task_profiler PINS (states EXEC/RELEASE_DEPS/
#     ACTIVATE_CB no .prof) + device_show_statistics + DAG (.dot via
#     CHAMELEON_PARSEC_DOT, so 8 runs pequenas -- o grapher nao pesa).
#
# Os modelos de desempenho do StarPU (dmda/dmdas) sao CALIBRADOS antes das
# medidas (CALIB_PASSES=2, modelo isolado em $RESULTS_DIR/.starpu, congelado
# p/ as runs cronometradas com STARPU_CALIBRATE=0).
#
# Normalmente disparado pelo front-end (scripts/pcad/collect_traces_tupi_pcad.sh),
# mas roda direto se voce ja esta no no:
#   bash scripts/pcad/collect_traces_tupi_node.sh
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
DOE_CSV="final/doe_gpu_traces_tupi_1000.csv"

if [[ -z "${SCRATCH:-}" ]]; then
    SCRATCH="$HOME/scratch"
    echo "!! \$SCRATCH nao definido; usando $SCRATCH (convencao do nixw)." >&2
    echo "   Prefira disparar via srun --jobid (herda o \$SCRATCH da alocacao)." >&2
fi

WORK_DIR="$SCRATCH/intro_hpc"
RESULTS_REL="data/gpu_traces_${NODE_TAG}_${JOBTAG}"
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

# DoE de traces para a tupi (8 runs: FP64, N=40000, b=1000, 1 rep).
cp "$DOE_CSV" "${RESULTS_REL}/design.csv"
export DESIGN_FILE="${RESULTS_REL}/design.csv"

# Rastro completo + contadores nativos de movimentacao de dados + DAG (.dot)
# do PaRSEC via CHAMELEON_PARSEC_DOT (so 8 runs pequenas, o grapher nao pesa).
export TRACE=1 TRACE_FULL=1 TRACE_DAG=1 TRACE_STATS=1
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
