#!/bin/bash

# Coleta MANUAL de traces, para rodar dentro de uma alocacao interativa (ex.:
# `sbatch --wrap 'sleep infinity'` no poti + ssh no no), sem passar por um job
# SLURM dedicado. Pega, para cada runtime x cada um dos seus dois escalonadores
# x cada algoritmo:
#   - StarPU (dmda, dmdas): trace FxT  (.../prof_file_*  -> StarVZ)
#   - PaRSEC (lfq, gd)    : trace .prof (PINS task_profiler) + DAG .dot
#
#   uso: bash scripts/collect_traces_manual.sh [starpu|parsec|both]
#
# Para o fluxo front-end -> no -> scratch -> home (igual aos jobs SLURM), use
# scripts/pcad/collect_traces_pcad.sh, que chama este script no scratch do no.
#
# Nao reimplementa a captura: apenas orquestra os dois dev shells (que sao
# *builds diferentes* do Chameleon, por isso nao da pra rodar os dois numa
# invocacao so) e delega o loop interno para scripts/gpu_doe_sweep.sh com os
# knobs de trace ligados. Mesma mecanica do slurm/full_trace.slurm, mas:
#   - sem o rsync pro $SCRATCH (roda no proprio repo);
#   - TRACE_DAG=1 (full_trace.slurm usa 0 para nao poluir o timeline grande);
#     aqui foi uma escolha consciente -- a DAG em n=19200 sai gigante.
#
# Env overrides: RESULTS_DIR, DESIGN_FILE, THREADS, GPUS, NIX_CMD.

set -euo pipefail

cd "$(dirname "$0")/.."   # raiz do repo

PHASE="${1:-both}"        # starpu | parsec | both

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS_DIR="${RESULTS_DIR:-data/manual_traces_${STAMP}}"
DESIGN_FILE="${DESIGN_FILE:-scripts/doe/doe_manual_traces.csv}"
THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"
# Na PCAD todo comando nix passa pelo shim nixw (ver CLAUDE.md). Localmente,
# rode com NIX_CMD="nix".
NIX_CMD="${NIX_CMD:-nixw nix}"

export RESULTS_DIR DESIGN_FILE THREADS GPUS
export RESULTS_FILE="${RESULTS_DIR}/results.csv"
export RUNS_DIR="${RESULTS_DIR}/runs"

# Knobs de trace consumidos por gpu_doe_sweep.sh:
export TRACE=1        # FxT (StarPU) / .prof (PaRSEC)
export TRACE_FULL=1   # PINS task_profiler -> .prof do PaRSEC com estados de runtime
export TRACE_DAG=1    # CHAMELEON_PARSEC_DOT -> DAG .dot (so PaRSEC)

if [[ ! -f "$DESIGN_FILE" ]]; then
    echo "DoE nao encontrado: $DESIGN_FILE" >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"
cp "$DESIGN_FILE" "$RESULTS_DIR/design.csv"

echo "==> Design : $DESIGN_FILE"
echo "==> Saida  : $RESULTS_DIR  (t=$THREADS g=$GPUS)"
echo "==> Fases  : $PHASE   (TRACE=$TRACE TRACE_FULL=$TRACE_FULL TRACE_DAG=$TRACE_DAG)"

# Roda uma fase: entra no dev shell <shell> e dispara o sweep filtrando <runtime>.
run_phase() {
    local shell="$1" runtime="$2"
    echo
    echo ">>> fase ${runtime}  (${NIX_CMD} develop .#${shell})"
    if ! $NIX_CMD develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" RESULTS_FILE="$RESULTS_FILE" RUNS_DIR="$RUNS_DIR" \
        DESIGN_FILE="$DESIGN_FILE" THREADS="$THREADS" GPUS="$GPUS" \
        TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" TRACE_DAG="$TRACE_DAG" \
        bash scripts/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! fase ${runtime} (.#${shell}) falhou; continuando" >&2
    fi
}

case "$PHASE" in
    starpu) run_phase starpu starpu ;;
    parsec) run_phase parsec parsec ;;
    both)   run_phase starpu starpu; run_phase parsec parsec ;;
    *) echo "uso: $0 [starpu|parsec|both]" >&2; exit 1 ;;
esac

echo
echo "==> Pronto. Traces por run em ${RUNS_DIR}/<idx>_<runtime>_<sched>_<kernel>_...:"
echo "      StarPU -> arquivos FxT (prof_file_*)"
echo "      PaRSEC -> cham_<kernel>-0.prof + cham_<kernel>.dot"
echo "    Resumo em ${RESULTS_FILE}"
