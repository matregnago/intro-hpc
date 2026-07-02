#!/usr/bin/env bash
#
# Build a self-contained StarPU-vs-PaRSEC comparison folder for validation.
# Runs every cross-runtime analysis script against ONE job (a job whose runs/
# holds both runtimes with the full StarVZ phase-1 parquet set) and routes all
# their figures/tables into a single separate folder via PLOTS_DIR.
#
# Produces, in $OUT:
#   st_compare_<algo>.{png,pdf}   side-by-side space-time (StarPU | PaRSEC)
#   bounds_makespan.{png,pdf}     makespan vs ABE/CPB/GPU bound, both runtimes
#   sched_concurrency.{png,pdf}   busy-worker concurrency profile, both runtimes
#   kiteration_compare_<algo>.*   k-iteration signature, all schedulers/runtimes
#   task_times_{mean,violin}.*    per-kernel duration, both runtimes
#   wait_time_dist.{png,pdf}      per-task DAG-readiness wait, both runtimes
#   transfers_d2h.{png,pdf}       H2D/D2H volume from log counters, both runtimes
#   *_summary.csv                 the backing tables
#
# Usage:  scripts/analysis/compare_all.sh [base_runs_dir] [out_dir]
#   base_runs_dir defaults to the GPU job (spotrf/sgeqrf n19200 b960).
#   out_dir       defaults to plots/compare.

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

BASE="${1:-data/manual_traces_20260621_172820/runs}"
OUT="${2:-plots/compare}"
mkdir -p "$OUT"
export PLOTS_DIR="$OUT"

echo "comparacao StarPU vs PaRSEC"
echo "  job:    $BASE"
echo "  saida:  $OUT"
echo

run() {  # run <label> <script> [args...]
  local label="$1"; shift
  echo "==> $label"
  if Rscript "$@"; then
    echo "    ok"
  else
    echo "    FALHOU (segue) -- $label" >&2
  fi
  echo
}

run "space-time (lado a lado)" "$HERE/plot_compare_st.r"      "$BASE"
run "bounds (makespan vs ABE/CPB/GPU)" "$HERE/plot_bounds.r"  "$BASE"
run "concurrency (escalonador)"  "$HERE/plot_scheduler_health.r" "$BASE"
run "k-iteration (por algoritmo)" "$HERE/plot_kiteration.r"   "$BASE"
run "task times (por kernel)"    "$HERE/plot_task_times.r"    "$BASE"
run "wait time (prontidao)"      "$HERE/plot_wait_time.r"     "$BASE"
run "transfers (H2D/D2H)"        "$HERE/plot_transfers.r"     "$BASE"

echo "pronto. figuras em: $OUT"
ls -1 "$OUT"
