#!/bin/bash

# DoE-driven GPU sweep for the StarPU-vs-PaRSEC scheduler study. Run *inside*
# the matching dev shell (.#starpu-cuda or .#parsec). Reads $DESIGN_FILE (a CSV
# produced by scripts/doe/doe_gpu_sched.r or doe_gpu_tile.r), runs the rows
# whose `runtime` matches the arg, and appends results to $RESULTS_FILE.
#
#   usage: gpu_doe_sweep.sh <starpu|parsec>
#
# Driven by scripts/run_gpu_sched.sh. DoE CSV columns (header skipped):
#   precision,algorithm,n,b,runtime,scheduler,rep
# `precision` in {FP32,FP64} picks the binary (chameleon_stesting/dtesting) and
# the kernel prefix (s/d); `algorithm` in {potrf,geqrf} -> e.g. spotrf, dgeqrf.
#
# Shares the validated GPU env of scripts/gpu_sweep.sh (see the
# local-gpu-compare memory): StarPU dmda/dmdas perf-model schedulers are
# calibrated once per (precision,algorithm) before the timed runs.

set -euo pipefail

RUNTIME_FILTER="${1:?usage: $0 <starpu|parsec>}"

RESULTS_DIR="${RESULTS_DIR:-results}"
RESULTS_FILE="${RESULTS_FILE:-${RESULTS_DIR}/results.csv}"
RUNS_DIR="${RUNS_DIR:-${RESULTS_DIR}/runs}"
DESIGN_FILE="${DESIGN_FILE:?DESIGN_FILE not set}"

THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"
# How many calibration passes to run per StarPU perf-model kernel before timing.
# One pass already yields hundreds of gemm/trsm/syrk samples (>> StarPU's
# CALIBRATE_MINIMUM=10), but potrf has only ~(n/b) tasks, so a couple of passes
# makes the model robust before we freeze it.
CALIB_PASSES="${CALIB_PASSES:-2}"

# Isolate the StarPU perf-model in a per-job dir so calibration is FRESH (not
# polluted by models from earlier experiments) and is captured alongside the
# results. dmda and dmdas share these history-based, codelet-keyed models.
export STARPU_HOME="${STARPU_HOME:-${RESULTS_DIR}/.starpu}"
mkdir -p "$STARPU_HOME"

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"
if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,scheduler,algorithm,precision,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"
fi

# Validated env (local-gpu-compare memory): one BLAS thread per worker, let
# StarPU bind freely (hwloc misdetects topology under WSL2).
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export STARPU_WORKERS_NOBIND=1
export STARPU_SILENT=1
export STARPU_WORKERS_GETBIND=0

bin_of()    { [[ "$1" == "FP32" ]] && echo chameleon_stesting || echo chameleon_dtesting; }
prefix_of() { [[ "$1" == "FP32" ]] && echo s || echo d; }

# StarPU perf-model schedulers (dmda/dmdas) need calibration. The models are
# keyed per codelet+tile-size and shared by dmda/dmdas, so we calibrate each
# unique (precision,algorithm,n,b) with CALIB_PASSES warm-up runs BEFORE the
# timed runs. These calibration runs are discarded (not appended to results).
# Afterwards the timed runs use STARPU_CALIBRATE=0 to FREEZE the model, so we
# measure steady-state scheduling, not model-fitting overhead.
calibrate_starpu() {
    export STARPU_CALIBRATE=1
    export STARPU_NCUDA="$GPUS"
    export STARPU_SCHED=dmda
    declare -A seen
    while IFS=',' read -r precision algorithm n b runtime scheduler rep; do
        [[ "$runtime" == "starpu" ]] || continue
        local key="${precision}_${algorithm}_${n}_${b}"
        [[ -n "${seen[$key]:-}" ]] && continue
        seen[$key]=1
        local bin kernel
        bin="$(bin_of "$precision")"
        kernel="$(prefix_of "$precision")${algorithm}"
        for pass in $(seq 1 "$CALIB_PASSES"); do
            echo "[$(date +%T)] calibrate starpu/$kernel n=$n b=$b (pass ${pass}/${CALIB_PASSES}) ..."
            "$bin" -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
                > /dev/null 2>&1 || true
        done
    done < <(tail -n +2 "$DESIGN_FILE")
}

[[ "$RUNTIME_FILTER" == "starpu" ]] && calibrate_starpu

idx=0
while IFS=',' read -r precision algorithm n b runtime scheduler rep; do
    idx=$((idx + 1))
    [[ "$runtime" == "$RUNTIME_FILTER" ]] || continue

    bin="$(bin_of "$precision")"
    kernel="$(prefix_of "$precision")${algorithm}"
    RUN_ID="$(printf '%04d' "$idx")_${runtime}_${scheduler}_${kernel}_n${n}_b${b}_rep${rep}"
    RUN_DIR="${RUNS_DIR}/${RUN_ID}"
    LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
    mkdir -p "$RUN_DIR"

    case "$runtime" in
        starpu)
            export STARPU_SCHED="$scheduler"
            export STARPU_NCUDA="$GPUS"
            export STARPU_CALIBRATE=0   # freeze the model: measure, don't fit
            unset PARSEC_MCA_mca_sched PARSEC_MCA_device_cuda_enabled || true
            ;;
        parsec)
            export PARSEC_MCA_mca_sched="$scheduler"
            export PARSEC_MCA_device_cuda_enabled=1
            unset STARPU_SCHED STARPU_NCUDA STARPU_CALIBRATE || true
            ;;
        *) echo "unknown runtime: $runtime" >&2; exit 1 ;;
    esac

    echo "[$(date +%T)] ${RUN_ID} (t=${THREADS} g=${GPUS})"
    set +e
    "$bin" -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
        > "$LOG_FILE" 2>&1
    status=$?
    set -e

    if line=$(grep -E '^[0-9]+;' "$LOG_FILE" | tail -n 1); [[ -n "${line:-}" ]]; then
        time=$(echo "$line"   | awk -F';' '{print $(NF-1)}')
        gflops=$(echo "$line" | awk -F';' '{print $NF}')
    else
        time="NA"; gflops="NA"
        echo "  !! no result line (exit ${status}); see ${LOG_FILE}" >&2
    fi

    echo "${runtime},${scheduler},${kernel},${precision},${n},${b},${THREADS},${GPUS},${rep},${time},${gflops}" \
        >> "$RESULTS_FILE"
done < <(tail -n +2 "$DESIGN_FILE")
