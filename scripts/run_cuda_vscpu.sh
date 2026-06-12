#!/bin/bash

# CPU-vs-GPU validation sweep for one runtime, run *inside* a dev shell.
#
#   usage: run_cuda_vscpu.sh <starpu|parsec> <cpu|gpu>
#
# Mirrors scripts/run.sh (same env knobs, same result-line parsing, same
# per-run log layout under runs/) but drives the cpu/gpu dimension and the
# small fixed sweep used to validate the CUDA backend. Unlike run.sh it does
# not read a DoE CSV; the sweep is defined by the env vars below so the same
# script covers all four (runtime x mode) phases of slurm/cuda_vs_cpu.slurm.

set -euo pipefail

RUNTIME="${1:?usage: $0 <starpu|parsec> <cpu|gpu>}"
MODE="${2:?usage: $0 <starpu|parsec> <cpu|gpu>}"

RESULTS_DIR="${RESULTS_DIR:-results}"
RESULTS_FILE="${RESULTS_DIR}/results.csv"
RUNS_DIR="${RESULTS_DIR}/runs"

# Sweep definition (overridable from the SLURM script / shell).
OPS="${OPS:-dpotrf dgeqrf}"
SWEEP_PAIRS="${SWEEP_PAIRS:-8000:1000 16000:1000 24000:1000}"
REPS="${REPS:-3}"

# Thread counts per mode. GPU uses fewer CPU workers because PaRSEC reserves
# one thread per GPU stream (see docs/parsec_cuda_results.md); 16/8 mirror the
# documented poti runs.
CPU_THREADS="${CPU_THREADS:-16}"
GPU_THREADS="${GPU_THREADS:-8}"

# One representative scheduler per runtime by default; pass space-separated
# lists to sweep more (e.g. STARPU_SCHEDS="lws ws").
STARPU_SCHEDS="${STARPU_SCHEDS:-lws}"
PARSEC_SCHEDS="${PARSEC_SCHEDS:-lfq}"

# VERBOSE=1 turns on PaRSEC's per-tile CUDA dispatch logging. Useful for the
# first validation run to confirm tiles land on the GPU; very noisy at large n.
VERBOSE="${VERBOSE:-0}"

case "$RUNTIME" in
    starpu) SCHEDS="$STARPU_SCHEDS" ;;
    parsec) SCHEDS="$PARSEC_SCHEDS" ;;
    *) echo "unknown runtime: $RUNTIME" >&2; exit 1 ;;
esac

case "$MODE" in
    cpu) threads="$CPU_THREADS"; gpus=0 ;;
    gpu) threads="$GPU_THREADS"; gpus=1 ;;
    *) echo "unknown mode: $MODE (expected cpu|gpu)" >&2; exit 1 ;;
esac

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"

if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,mode,scheduler,algorithm,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export STARPU_WORKERS_NOBIND=1
export STARPU_SILENT=1

idx=0
for scheduler in $SCHEDS; do
    # Runtime-specific environment. The cpu/gpu shells differ in the binary
    # (CUDA-enabled or not); here we only toggle the runtime's scheduler and,
    # for PaRSEC, whether the CUDA device is enabled at runtime.
    case "$RUNTIME" in
        starpu)
            export STARPU_SCHED="$scheduler"
            unset PARSEC_MCA_mca_sched PARSEC_MCA_device_cuda_enabled \
                  PARSEC_MCA_device_show_capabilities \
                  PARSEC_MCA_device_cuda_verbose PARSEC_MCA_device_verbose
            ;;
        parsec)
            export PARSEC_MCA_mca_sched="$scheduler"
            unset STARPU_SCHED
            if [[ "$MODE" == "gpu" ]]; then
                export PARSEC_MCA_device_cuda_enabled=1
                export PARSEC_MCA_device_show_capabilities=1
                if [[ "$VERBOSE" == "1" ]]; then
                    export PARSEC_MCA_device_cuda_verbose=20
                    export PARSEC_MCA_device_verbose=20
                else
                    unset PARSEC_MCA_device_cuda_verbose PARSEC_MCA_device_verbose
                fi
            else
                export PARSEC_MCA_device_cuda_enabled=0
                unset PARSEC_MCA_device_show_capabilities \
                      PARSEC_MCA_device_cuda_verbose PARSEC_MCA_device_verbose
            fi
            ;;
    esac

    for algorithm in $OPS; do
        for pair in $SWEEP_PAIRS; do
            n="${pair%%:*}"
            b="${pair##*:}"
            for rep in $(seq 1 "$REPS"); do
                idx=$((idx + 1))
                RUN_ID="$(printf '%04d' "$idx")_${RUNTIME}_${MODE}_${scheduler}_${algorithm}_n${n}_b${b}_rep${rep}"
                RUN_DIR="${RUNS_DIR}/${RUN_ID}"
                LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
                mkdir -p "$RUN_DIR"

                echo "[$(date +%T)] ${RUN_ID} (t=${threads} g=${gpus})"

                set +e
                chameleon_stesting \
                    -o "$algorithm" \
                    -n "$n" \
                    -b "$b" \
                    -t "$threads" \
                    -g "$gpus" \
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

                echo "${RUNTIME},${MODE},${scheduler},${algorithm},${n},${b},${threads},${gpus},${rep},${time},${gflops}" \
                    >> "$RESULTS_FILE"
            done
        done
    done
done
