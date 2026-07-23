#!/bin/bash

set -euo pipefail

RUNTIME=$1
RESULTS_DIR="${RESULTS_DIR:-results}"
DESIGN_FILE="${DESIGN_FILE:?DESIGN_FILE not set}"
RESULTS_FILE="${RESULTS_DIR}/results.csv"
RUNS_DIR="${RESULTS_DIR}/runs"

THREADS="${THREADS:-12}"
GPUS="${GPUS:-1}"
TRACE="${TRACE:-0}"
TRACE_FULL="${TRACE_FULL:-0}"
TRACE_DAG="${TRACE_DAG:-$TRACE_FULL}"
TRACE_STATS="${TRACE_STATS:-$TRACE}"
CALIB_PASSES="${CALIB_PASSES:-2}"

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"

if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,scheduler,algorithm,precision,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Fresh StarPU perf-models, stored alongside the results.
export STARPU_HOME="${RESULTS_DIR}/.starpu"
mkdir -p "$STARPU_HOME"

# CUDA streams per GPU: match PaRSEC's default of 4 kernel streams.
export STARPU_NWORKER_PER_CUDA="${STARPU_NWORKER_PER_CUDA:-4}"

# nixglhos injects the host's GPU drivers.
chameleon() { nixglhost chameleon_dtesting "$@"; }

calibrate_starpu() {
    export STARPU_CALIBRATE=1
    export STARPU_NCUDA="$GPUS"
    export STARPU_SCHED=dmda

    declare -A seen
    local n_combos=0
    while IFS=',' read -r precision algorithm n b runtime scheduler rep; do
        [[ "$runtime" == "starpu" ]] || continue
        local key="${precision}_${algorithm}_${n}_${b}"
        [[ -n "${seen[$key]:-}" ]] && continue
        seen[$key]=1
        n_combos=$((n_combos + 1))

        local kernel="d${algorithm}"
        for pass in $(seq 1 "$CALIB_PASSES"); do
            echo "[$(date +%T)] calibrate starpu/${kernel} n=${n} b=${b} (pass ${pass}/${CALIB_PASSES})"
            chameleon -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
                > /dev/null 2>&1 || true
        done
    done < <(tail -n +2 "$DESIGN_FILE")
}

if [[ "$RUNTIME" == "starpu" && "$CALIB_PASSES" -gt 0 ]]; then
    calibrate_starpu
fi

EXTRA_FLAGS=()
case "$RUNTIME" in
    starpu)
        export STARPU_NCUDA="$GPUS"
        export STARPU_CALIBRATE=0
        export STARPU_FXT_TRACE="$TRACE"
        if [[ "$TRACE_STATS" == "1" ]]; then
            export STARPU_BUS_STATS=1 STARPU_WORKER_STATS=1
        fi
        if [[ "$TRACE" == "1" ]]; then
            EXTRA_FLAGS+=(--trace)
        fi
        ;;
    parsec)
        export PARSEC_MCA_device_cuda_enabled=1
        if [[ "$TRACE_STATS" == "1" ]]; then
            export PARSEC_MCA_device_show_statistics=1
        fi
        if [[ "$TRACE" == "1" && "$TRACE_FULL" == "1" ]]; then
            export PARSEC_MCA_mca_pins=task_profiler
        fi
        ;;
    *)
        exit 1
        ;;
esac

idx=0
while IFS=',' read -r precision algorithm n b runtime scheduler rep; do
    idx=$((idx + 1))
    [[ "$runtime" == "$RUNTIME" ]] || continue

    kernel="d${algorithm}"
    RUN_ID="$(printf '%04d' "$idx")_${runtime}_${scheduler}_${kernel}_n${n}_b${b}_rep${rep}"
    RUN_DIR="${RUNS_DIR}/${RUN_ID}"
    LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
    mkdir -p "$RUN_DIR"

    echo "[$(date +%T)] ${RUN_ID} (t=${THREADS} g=${GPUS} trace=${TRACE})"

    if [[ "$RUNTIME" == "starpu" ]]; then
        export STARPU_SCHED="$scheduler"
        if [[ "$TRACE" == "1" ]]; then
            export STARPU_FXT_PREFIX="${RUN_DIR}/"
        fi
    else
        export PARSEC_MCA_mca_sched="$scheduler"
        if [[ "$TRACE" == "1" ]]; then
            export PARSEC_MCA_profile_filename="${RUN_DIR}/cham_${kernel}"
            if [[ "$TRACE_DAG" == "1" ]]; then
                export CHAMELEON_PARSEC_DOT="${RUN_DIR}/cham_${kernel}"
            fi
        fi
    fi

    set +e
    chameleon -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
        "${EXTRA_FLAGS[@]}" > "$LOG_FILE" 2>&1
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

    echo "${runtime},${scheduler},${kernel},${precision},${n},${b},${THREADS},${GPUS},${rep},${time},${gflops}" \
        >> "$RESULTS_FILE"

done < <(tail -n +2 "$DESIGN_FILE")
