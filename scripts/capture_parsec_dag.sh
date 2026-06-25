#!/bin/bash

set -euo pipefail

DESIGN_FILE="${DESIGN_FILE:-scripts/doe/doe_parsec_dag.csv}"
RESULTS_DIR="${RESULTS_DIR:-data/parsec_dag}"
RUNS_DIR="$RESULTS_DIR/runs"
THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"

mkdir -p "$RUNS_DIR"

# 1 thread BLAS por worker; o resto dos nucleos fica para o PaRSEC + a GPU.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PARSEC_MCA_device_cuda_enabled="$GPUS"
export PARSEC_MCA_mca_pins=task_profiler

while IFS=',' read -r algorithm n b scheduler; do
    run_dir="$RUNS_DIR/${algorithm}_${scheduler}_n${n}_b${b}"
    mkdir -p "$run_dir"

    export PARSEC_MCA_mca_sched="$scheduler"
    export PARSEC_MCA_profile_filename="$run_dir/cham_${algorithm}"
    export CHAMELEON_PARSEC_DOT="$run_dir/cham_${algorithm}"

    echo "[$(date +%T)] ${algorithm} sched=${scheduler} n=${n} b=${b}"
    chameleon_stesting -o "$algorithm" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" \
        --nowarmup > "$run_dir/run.log" 2>&1 \
        || echo "  !! falhou; ver $run_dir/run.log" >&2
done < <(tail -n +2 "$DESIGN_FILE")

echo "==> DAGs (.dot) e profiles (.prof) em $RUNS_DIR/"
