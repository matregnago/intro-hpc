#!/bin/bash

# Roda NO NO de computacao (poti), dentro da sua alocacao interativa
# (`sbatch --wrap 'sleep infinity'`). E o equivalente manual do
# slurm/full_trace.slurm: rsync do repo -> $SCRATCH, roda as duas fases de
# trace la (I/O rapido, nix-portable ja instalado em $SCRATCH/nix-data), e
# copia os resultados de volta pro $HOME via trap EXIT (copyback). A captura em
# si fica em scripts/collect_traces_manual.sh -- aqui so cuidamos do scratch.
#
# Normalmente disparado pelo front-end (scripts/pcad/collect_traces_pcad.sh),
# mas roda direto se voce ja esta no no:
#   bash scripts/pcad/collect_traces_node.sh
#
# Env: SUBMIT_DIR (repo no home, default ~/intro_hpc), SCRATCH, JOBTAG,
#      THREADS, GPUS, PHASE (starpu|parsec|both).

set -euo pipefail

SUBMIT_DIR="${SUBMIT_DIR:-$HOME/intro_hpc}"
JOBTAG="${JOBTAG:-$(date +%Y%m%d_%H%M%S)}"
PHASE="${PHASE:-both}"

# $SCRATCH normalmente vem da alocacao SLURM (igual aos jobs). Se rodou via ssh
# cru e ele nao veio, usamos a mesma convencao do nixw ($HOME/scratch) -- mas
# avisamos, porque um scratch diferente do da alocacao faz o nixw rebaixar e
# rebaixar/rebaixar o store nix-portable.
if [[ -z "${SCRATCH:-}" ]]; then
    SCRATCH="$HOME/scratch"
    echo "!! \$SCRATCH nao definido; usando $SCRATCH (convencao do nixw)." >&2
    echo "   Prefira disparar via srun --jobid (herda o \$SCRATCH da alocacao)." >&2
fi

WORK_DIR="$SCRATCH/intro_hpc"
RESULTS_REL="data/manual_traces_${JOBTAG}"
LOCAL_RESULTS_DIR="$WORK_DIR/$RESULTS_REL"
FINAL_RESULTS_DIR="$SUBMIT_DIR/$RESULTS_REL"

echo "==============================================================="
echo "Node       : $(hostname)"
echo "SubmitDir  : $SUBMIT_DIR"
echo "WorkDir    : $WORK_DIR   (\$SCRATCH=$SCRATCH)"
echo "JobTag     : $JOBTAG"
echo "==============================================================="

[[ -f "$SUBMIT_DIR/flake.nix" ]] || { echo "repo nao encontrado em $SUBMIT_DIR" >&2; exit 1; }

mkdir -p "$WORK_DIR" "$FINAL_RESULTS_DIR"
# Espelha o repo no scratch. Exclui os dados pesados (a captura escreve um data/
# novo no scratch) e o .git/.direnv. --delete mantem o WORK_DIR limpo, mas
# 'data/' fica excluido do --delete, entao traces de runs anteriores no scratch
# nao sao apagados.
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
export THREADS="${THREADS:-16}" GPUS="${GPUS:-1}"

bash scripts/collect_traces_manual.sh "$PHASE"

echo "==> Traces em $FINAL_RESULTS_DIR  (no \$HOME, compartilhado com o front-end)"
