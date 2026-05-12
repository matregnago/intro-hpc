#!/bin/bash

set -euo pipefail

RUNTIME_FILTER="${1:?usage: $0 <starpu|parsec>}"
RESULTS_DIR="${RESULTS_DIR:-results}"
DESIGN_FILE="${DESIGN_FILE:-doe_n_size.csv}"
RESULTS_FILE="${RESULTS_DIR}/results.csv"
RUNS_DIR="${RESULTS_DIR}/runs"

TRACE="${TRACE:-0}"

threads=24

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"

if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,scheduler,algorithm,n,b,threads,rep,time,gflops" > "$RESULTS_FILE"
fi

export OPENBLAS_NUM_THREADS=1
export STARPU_WORKERS_NOBIND=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export STARPU_SILENT=1

idx=0
while IFS=',' read -r -a fields; do
    idx=$((idx + 1))
    algorithm="${fields[0]}"
    n="${fields[1]}"
    b="${fields[2]}"
    runtime="${fields[3]}"
    scheduler="${fields[4]}"
    rep="${fields[-1]}"
    if [[ "$runtime" != "$RUNTIME_FILTER" ]]; then
        continue
    fi

    rep_tag="${rep#.}"
    RUN_ID="$(printf '%04d' "$idx")_${runtime}_${scheduler}_${algorithm}_n${n}_b${b}_rep${rep_tag}"
    RUN_DIR="${RUNS_DIR}/${RUN_ID}"
    LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
    TRACE_DIR="$RUN_DIR"
    mkdir -p "$RUN_DIR"

    echo "[$(date +%T)] ${RUN_ID}"

    EXTRA_FLAGS=()
    case "$runtime" in
        starpu)
            export STARPU_SCHED="$scheduler"
            export STARPU_NCPU="$threads"
            unset PARSEC_MCA_mca_sched PARSEC_MCA_profile_filename
            if [[ "$TRACE" == "1" ]]; then
                export STARPU_FXT_TRACE=1
                export STARPU_FXT_PREFIX="${TRACE_DIR}/"
                EXTRA_FLAGS+=(--trace)
            else
                export STARPU_FXT_TRACE=0
                unset STARPU_FXT_PREFIX
            fi
            ;;
        parsec)
            export PARSEC_MCA_mca_sched="$scheduler"
            unset STARPU_SCHED STARPU_NCPU STARPU_FXT_TRACE STARPU_FXT_PREFIX
            if [[ "$TRACE" == "1" ]]; then
                export PARSEC_MCA_profile_filename="${TRACE_DIR}/cham_${algorithm}"
            else
                unset PARSEC_MCA_profile_filename
            fi
            ;;
        *)
            echo "unknown runtime: $runtime" >&2
            exit 1
            ;;
    esac

    set +e
    chameleon_dtesting \
        -o "$algorithm" \
        -n "$n" \
        -b "$b" \
        -t "$threads" \
        "${EXTRA_FLAGS[@]}" \
        > "$LOG_FILE" 2>&1
    status=$?
    set -e

    if line=$(grep -E '^[0-9]+;' "$LOG_FILE" | tail -n 1); [[ -n "${line:-}" ]]; then
        time=$(echo "$line"   | awk -F';' '{print $(NF-1)}')
        gflops=$(echo "$line" | awk -F';' '{print $NF}')
    else
        time="NA"
        gflops="NA"
    fi

    echo "${runtime},${scheduler},${algorithm},${n},${b},${threads},${rep},${time},${gflops}" \
        >> "$RESULTS_FILE"

done < <(tail -n +2 "$DESIGN_FILE")
