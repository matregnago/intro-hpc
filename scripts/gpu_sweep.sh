#!/bin/bash

# Per-runtime GPU sweep, run *inside* the matching dev shell
# (.#starpu-cuda or .#parsec). Sweeps {dpotrf,dgeqrf} x N x reps for one
# runtime with the GPU enabled (-g 1) and appends rows to $RESULTS_FILE.
#
#   usage: gpu_sweep.sh <starpu|parsec>
#
# Driven by scripts/run_local_gpu_compare.sh, which enters each dev shell and
# defines the sweep via the env vars below. Mirrors the env knobs of
# scripts/run.sh / run_cuda_vscpu.sh but always double precision
# (chameleon_dtesting) and GPU-only, since the goal is a StarPU-vs-PaRSEC
# comparison with CUDA on, not a CPU baseline.

set -euo pipefail

RUNTIME="${1:?usage: $0 <starpu|parsec>}"

RESULTS_DIR="${RESULTS_DIR:-results}"
RESULTS_FILE="${RESULTS_FILE:-${RESULTS_DIR}/results.csv}"
RUNS_DIR="${RUNS_DIR:-${RESULTS_DIR}/runs}"

# Sweep definition (overridable from the orchestrator / shell).
# Single precision by default (spotrf/sgeqrf via chameleon_stesting): FP32 is
# where this consumer GPU actually beats the CPU; FP64 (dpotrf/dgeqrf via
# chameleon_dtesting) is hardware-bound below the 16-thread CPU on an RTX 4060
# Ti. To run double instead: TESTING_BIN=chameleon_dtesting OPS="dpotrf dgeqrf".
TESTING_BIN="${TESTING_BIN:-chameleon_stesting}"
OPS="${OPS:-spotrf sgeqrf}"
NS="${NS:-9600 19200 28800 38400}"
B="${B:-960}"
THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"
REPS="${REPS:-3}"

# Best-effort GPU scheduler per runtime. StarPU's lws explicitly warns it is
# not GPU-aware; dmda (perf-model based) roughly doubles GPU throughput here and
# does NOT hang with --nowarmup, so we calibrate it (model persists in
# ~/.starpu, keyed per kernel+tile-size, so b is fixed across the sweep).
# PaRSEC has no perf-model scheduler; lfq is its standard ready-task queue and
# GPU placement is decided by PaRSEC's heterogeneous device heuristic regardless.
STARPU_SCHED_NAME="${STARPU_SCHED_NAME:-dmda}"
PARSEC_SCHED_NAME="${PARSEC_SCHED_NAME:-lfq}"

case "$RUNTIME" in
    starpu) scheduler="$STARPU_SCHED_NAME" ;;
    parsec) scheduler="$PARSEC_SCHED_NAME" ;;
    *) echo "unknown runtime: $RUNTIME" >&2; exit 1 ;;
esac

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"

# Header is written once by the orchestrator; create it here too so the script
# is usable standalone.
if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,scheduler,algorithm,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"
fi

# One BLAS thread per worker: the runtime owns the parallelism.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export STARPU_WORKERS_NOBIND=1
export STARPU_SILENT=1
# hwloc under WSL2 misdetects the topology (reports 10 cores/20 threads) and
# StarPU oversubscribes PUs; let StarPU bind freely instead.
export STARPU_WORKERS_GETBIND=0

# Runtime-specific environment (the dev shell already provides the right
# CUDA-enabled chameleon_dtesting binary and the WSL libcuda driver path).
case "$RUNTIME" in
    starpu)
        export STARPU_SCHED="$scheduler"
        export STARPU_NCUDA="$GPUS"
        # Perf-model schedulers (dmda/dmdas) need calibration; the model
        # persists in ~/.starpu and is reused across the sweep.
        export STARPU_CALIBRATE="${STARPU_CALIBRATE:-1}"
        unset PARSEC_MCA_mca_sched PARSEC_MCA_device_cuda_enabled || true
        ;;
    parsec)
        export PARSEC_MCA_mca_sched="$scheduler"
        export PARSEC_MCA_device_cuda_enabled=1
        unset STARPU_SCHED STARPU_NCUDA || true
        ;;
esac

# StarPU dmda warmup: calibrate each op's perf models once at the smallest N
# (models are per kernel + tile-size, and b is fixed, so this primes every
# kernel for the whole sweep). Discarded — not recorded.
if [[ "$RUNTIME" == "starpu" ]]; then
    warm_n=$(echo "$NS" | tr ' ' '\n' | sort -n | head -1)
    for algorithm in $OPS; do
        echo "[$(date +%T)] calibrating starpu/$algorithm at n=$warm_n ..."
        "$TESTING_BIN" -o "$algorithm" -n "$warm_n" -b "$B" -t "$THREADS" \
            -g "$GPUS" --nowarmup > /dev/null 2>&1 || true
    done
fi

idx=0
for algorithm in $OPS; do
    for n in $NS; do
        for rep in $(seq 1 "$REPS"); do
            idx=$((idx + 1))
            RUN_ID="$(printf '%04d' "$idx")_${RUNTIME}_${scheduler}_${algorithm}_n${n}_b${B}_rep${rep}"
            RUN_DIR="${RUNS_DIR}/${RUN_ID}"
            LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
            mkdir -p "$RUN_DIR"

            echo "[$(date +%T)] ${RUN_ID} (t=${THREADS} g=${GPUS})"

            set +e
            "$TESTING_BIN" \
                -o "$algorithm" \
                -n "$n" \
                -b "$B" \
                -t "$THREADS" \
                -g "$GPUS" \
                --nowarmup \
                > "$LOG_FILE" 2>&1
            status=$?
            set -e

            if line=$(grep -E '^[0-9]+;' "$LOG_FILE" | tail -n 1); [[ -n "${line:-}" ]]; then
                time=$(echo "$line"   | awk -F';' '{print $(NF-1)}')
                gflops=$(echo "$line" | awk -F';' '{print $NF}')
            else
                time="NA"
                gflops="NA"
                echo "  !! no result line (exit ${status}); see ${LOG_FILE}" >&2
            fi

            echo "${RUNTIME},${scheduler},${algorithm},${n},${B},${THREADS},${GPUS},${rep},${time},${gflops}" \
                >> "$RESULTS_FILE"
        done
    done
done
