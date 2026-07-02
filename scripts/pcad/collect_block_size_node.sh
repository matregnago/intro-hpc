#!/bin/bash

# Roda NO NO de computacao (poti ou tupi), dentro da sua alocacao interativa
# (`sbatch -p <poti|tupi> --wrap 'sleep infinity'`). E o equivalente manual do
# slurm/gpu_traces_fair.slurm mas para a fase de SELECAO DE TILE (E1):
# rsync do repo -> $SCRATCH, roda as duas fases (starpu/parsec) do
# gpu_doe_sweep.sh la (I/O rapido), e copia os resultados de volta pro $HOME
# via trap EXIT (copyback). Sem traces (so time/gflops) -- TRACE=0.
#
# O no e detectado via hostname: poti -> THREADS=19 + doe_gpu_tile_poti.csv,
# tupi -> THREADS=23 + doe_gpu_tile_tupi.csv.
#
# Normalmente disparado pelo front-end (scripts/pcad/collect_block_size_pcad.sh),
# mas roda direto se voce ja esta no no:
#   bash scripts/pcad/collect_block_size_node.sh
#
# Env: SUBMIT_DIR (repo no home, default ~/intro_hpc), SCRATCH, JOBTAG,
#      THREADS (override), GPUS (1).

set -euo pipefail

SUBMIT_DIR="${SUBMIT_DIR:-$HOME/intro_hpc}"
JOBTAG="${JOBTAG:-$(date +%Y%m%d_%H%M%S)}"

# Detecta o no pelo hostname e seleciona CSV + THREADS. poti=24 cores (19 p/
# workers, 1 GPU driver + 4 PaRSEC streams), tupi=24 cores (23 p/ workers).
HOST="$(hostname)"
case "$HOST" in
    *poti*) NODE_TAG="poti"; DEFAULT_THREADS=19 ;;
    *tupi*) NODE_TAG="tupi"; DEFAULT_THREADS=23 ;;
    *)    echo "!! no desconhecido: $HOST (esperado poti ou tupi)" >&2; exit 1 ;;
esac
DOE_CSV="scripts/doe/doe_gpu_tile_${NODE_TAG}.csv"

# $SCRATCH normalmente vem da alocacao SLURM. Se rodou via ssh cru e ele nao
# veio, usamos a mesma convencao do nixw ($HOME/scratch).
if [[ -z "${SCRATCH:-}" ]]; then
    SCRATCH="$HOME/scratch"
    echo "!! \$SCRATCH nao definido; usando $SCRATCH (convencao do nixw)." >&2
    echo "   Prefira disparar via srun --jobid (herda o \$SCRATCH da alocacao)." >&2
fi

WORK_DIR="$SCRATCH/intro_hpc"
RESULTS_REL="data/gpu_tile_${NODE_TAG}_${JOBTAG}"
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

# DoE de tile para este no (FP64, 5 reps, 120 runs).
cp "$DOE_CSV" "${RESULTS_REL}/design.csv"
export DESIGN_FILE="${RESULTS_REL}/design.csv"

# Sem traces nesta fase -- so medimos time/gflops.
export TRACE=0

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> phase: ${runtime}  (nixw nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" TRACE="$TRACE" \
        bash scripts/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! phase ${runtime} (.#${shell}) failed; continuing" >&2
    fi
}

run_phase starpu starpu
run_phase parsec parsec

echo "==> Resultados em $FINAL_RESULTS_DIR  (no \$HOME, compartilhado com o front-end)"
