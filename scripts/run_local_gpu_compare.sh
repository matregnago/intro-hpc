#!/bin/bash

# Local StarPU-vs-PaRSEC GPU comparison (mini-simulation to validate the
# environment + scripts before running the real sweep on PCAD).
#
#   usage: bash scripts/run_local_gpu_compare.sh
#
# Runs the same {dpotrf,dgeqrf} x N sweep for both runtimes with the GPU
# enabled, entering each runtime's CUDA dev shell once, then renders the
# comparison plots. Tunable via env (defaults in scripts/gpu_sweep.sh):
#   OPS, NS, B, THREADS, GPUS, REPS
#
# Results land in data/local_gpu_compare_<timestamp>/ :
#   results.csv            one row per run
#   runs/<id>/<id>.log     per-run chameleon log
#   plots/                 generated figures + summary
#
# Both runtimes write to the same results.csv (they append). The plot step
# runs in the default dev shell (which has R + ggplot).

set -euo pipefail

cd "$(dirname "$0")/.."   # repo root

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS_DIR="${RESULTS_DIR:-data/local_gpu_compare_${STAMP}}"
export RESULTS_DIR
export RESULTS_FILE="${RESULTS_DIR}/results.csv"
export RUNS_DIR="${RESULTS_DIR}/runs"

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"
echo "runtime,scheduler,algorithm,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"

# Forward the sweep knobs into each dev shell.
PASS_ENV=(RESULTS_DIR RESULTS_FILE RUNS_DIR TESTING_BIN OPS NS B THREADS GPUS REPS
          STARPU_SCHED_NAME PARSEC_SCHED_NAME)
env_args=()
for v in "${PASS_ENV[@]}"; do
    [[ -n "${!v:-}" ]] && env_args+=("$v=${!v}")
done

echo "==> Results dir: $RESULTS_DIR"
echo "==> Sweep: BIN='${TESTING_BIN:-chameleon_dtesting}' OPS='${OPS:-dpotrf dgeqrf}' NS='${NS:-9600 19200 28800 38400}' B='${B:-960}' THREADS='${THREADS:-16}' REPS='${REPS:-3}'"

echo
echo "===== StarPU (.#starpu) ====="
nix develop .#starpu --impure --command env "${env_args[@]}" bash scripts/gpu_sweep.sh starpu

echo
echo "===== PaRSEC (.#parsec) ====="
nix develop .#parsec --impure --command env "${env_args[@]}" bash scripts/gpu_sweep.sh parsec

echo
echo "===== Plots (default shell) ====="
nix develop --impure --command Rscript scripts/analysis/plot_local_gpu_compare.r "$RESULTS_DIR"

echo
echo "==> Done. See ${RESULTS_DIR}/plots/"
