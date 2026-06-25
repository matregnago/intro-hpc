#!/bin/bash

# Orchestrator for the StarPU-vs-PaRSEC GPU *scheduler-impact* study (E2), the
# single-node GPU analogue of Schnorr's WAMCA-2025 nmad/openmpi experiment.
#
#   usage: bash scripts/run_gpu_sched.sh
#
# Reads a DoE CSV ($DESIGN_FILE, default scripts/doe/doe_gpu_sched.csv), runs
# every row of it once per runtime inside that runtime's CUDA dev shell, then
# renders the comparison plots. Tunable via env: DESIGN_FILE, THREADS, GPUS.
#
# Output: data/gpu_sched_<timestamp>/{results.csv,design.csv,runs/,plots/}.
#
# NOTE: uses plain `nix`. On PCAD prefix with the nixw shim (see CLAUDE.md) or
# adapt into a slurm/ job like the existing ones.

set -euo pipefail

cd "$(dirname "$0")/.."   # repo root

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS_DIR="${RESULTS_DIR:-data/gpu_sched_${STAMP}}"
export RESULTS_DIR
export RESULTS_FILE="${RESULTS_DIR}/results.csv"
export RUNS_DIR="${RESULTS_DIR}/runs"
export DESIGN_FILE="${DESIGN_FILE:-scripts/doe/doe_gpu_sched.csv}"

if [[ ! -f "$DESIGN_FILE" ]]; then
    echo "missing $DESIGN_FILE -- generate it first:" >&2
    echo "  nix develop --impure --command Rscript scripts/doe/doe_gpu_sched.r" >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"
cp "$DESIGN_FILE" "$RESULTS_DIR/design.csv"
echo "runtime,scheduler,algorithm,precision,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"

PASS_ENV=(RESULTS_DIR RESULTS_FILE RUNS_DIR DESIGN_FILE THREADS GPUS TRACE TRACE_FULL)
env_args=()
for v in "${PASS_ENV[@]}"; do
    [[ -n "${!v:-}" ]] && env_args+=("$v=${!v}")
done

echo "==> Design : $DESIGN_FILE"
echo "==> Results: $RESULTS_DIR"

echo
echo "===== StarPU (.#starpu-cuda) ====="
nix develop .#starpu-cuda --impure --command env "${env_args[@]}" bash scripts/gpu_doe_sweep.sh starpu

echo
echo "===== PaRSEC (.#parsec) ====="
nix develop .#parsec --impure --command env "${env_args[@]}" bash scripts/gpu_doe_sweep.sh parsec

echo
echo "===== Plots (default shell) ====="
nix develop --impure --command Rscript scripts/analysis/plot_gpu_sched.r "$RESULTS_DIR"

echo
echo "==> Done. See ${RESULTS_DIR}/plots/"
