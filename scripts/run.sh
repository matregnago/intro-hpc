#!/bin/bash

set -euo pipefail

RUNTIME=$1
RESULTS_DIR="${RESULTS_DIR:-results}"
DESIGN_FILE="${DESIGN_FILE:?DESIGN_FILE not set}"
RESULTS_FILE="${RESULTS_DIR}/results.csv"
RUNS_DIR="${RESULTS_DIR}/runs"

THREADS="${THREADS:-12}"
GPUS="${GPUS:-1}"
# TRACE=1 captures StarPU FxT / PaRSEC profile files next to each run's log.
TRACE="${TRACE:-0}"
# TRACE_FULL=1 (needs TRACE=1) adds PaRSEC's task_profiler PINS module, so the
# .prof also carries runtime/scheduler states (EXEC/RELEASE_DEPS/ACTIVATE_CB).
# StarPU's FxT is already complete under TRACE=1, no extra knob needed there.
TRACE_FULL="${TRACE_FULL:-0}"
# TRACE_DAG=1 (PaRSEC only) emits the task-graph DAG (.dot) via
# CHAMELEON_PARSEC_DOT. Decoupled from TRACE_FULL: at large N the DAG is a huge
# unrenderable hairball and the grapher pollutes the timed trace, so it's meant
# for a small dedicated run. Defaults to TRACE_FULL so simple callers (just
# TRACE_FULL=1) keep getting the DAG unless they opt out explicitly.
TRACE_DAG="${TRACE_DAG:-$TRACE_FULL}"
# TRACE_STATS=1 prints native data-movement counters into each run's log
# (StarPU STARPU_BUS_STATS/WORKER_STATS; PaRSEC device_show_statistics) --
# independent of the trace itself, no measurable overhead, so default it on
# whenever TRACE is on.
TRACE_STATS="${TRACE_STATS:-$TRACE}"
# StarPU perf-model calibration passes per unique (precision,algorithm,n,b)
# before the timed runs (dmda/dmdas are model-based schedulers; a cold model on
# the first timed rep would bias the measurement). Set to 0 to skip.
CALIB_PASSES="${CALIB_PASSES:-2}"

mkdir -p "$RESULTS_DIR" "$RUNS_DIR"

if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "runtime,scheduler,algorithm,precision,n,b,threads,gpus,rep,time,gflops" > "$RESULTS_FILE"
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Isolate the StarPU perf-model in this job's results dir, so calibration is
# fresh (not polluted by models from earlier experiments) and travels with the
# results. dmda and dmdas share these history-based, codelet-keyed models.
export STARPU_HOME="${RESULTS_DIR}/.starpu"
mkdir -p "$STARPU_HOME"

# Workers (CUDA streams) per GPU. StarPU defaults to 1 execution stream per
# device; PaRSEC defaults to 4 (device_cuda_max_streams=6 -> 2 transfer + 4
# kernel), so match it at 4 for a paired CPU+GPU comparison. No-op on CPU-only
# runs; ignored by PaRSEC.
export STARPU_NWORKER_PER_CUDA="${STARPU_NWORKER_PER_CUDA:-4}"

# precision (FP32/FP64) picks the test binary and the kernel prefix for `-o`.
bin_of()    { [[ "$1" == "FP32" ]] && echo chameleon_stesting || echo chameleon_dtesting; }
prefix_of() { [[ "$1" == "FP32" ]] && echo s || echo d; }

# StarPU perf-model schedulers (dmda/dmdas) need calibration. Models are keyed
# per codelet+tile-size and SHARED by dmda/dmdas, so we calibrate once per
# unique (precision,algorithm,n,b) -- always under dmda -- with CALIB_PASSES
# warm-up runs BEFORE the timed loop, discarded (not appended to results).
# Afterwards the timed loop freezes the model (STARPU_CALIBRATE=0).
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

        local bin kernel
        bin="$(bin_of "$precision")"
        kernel="$(prefix_of "$precision")${algorithm}"
        for pass in $(seq 1 "$CALIB_PASSES"); do
            echo "[$(date +%T)] calibrate starpu/${kernel} n=${n} b=${b} (pass ${pass}/${CALIB_PASSES})"
            "$bin" -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
                > /dev/null 2>&1 || true
        done
    done < <(tail -n +2 "$DESIGN_FILE")
    echo "[calib] concluido: ${n_combos} combinacoes (precision,algorithm,n,b)"
}

if [[ "$RUNTIME" == "starpu" && "$CALIB_PASSES" -gt 0 ]]; then
    calibrate_starpu
fi

idx=0
while IFS=',' read -r precision algorithm n b runtime scheduler rep; do
    idx=$((idx + 1))
    [[ "$runtime" == "$RUNTIME" ]] || continue

    bin="$(bin_of "$precision")"
    kernel="$(prefix_of "$precision")${algorithm}"
    rep_tag="${rep#.}"
    RUN_ID="$(printf '%04d' "$idx")_${runtime}_${scheduler}_${kernel}_n${n}_b${b}_rep${rep_tag}"
    RUN_DIR="${RUNS_DIR}/${RUN_ID}"
    LOG_FILE="${RUN_DIR}/${RUN_ID}.log"
    mkdir -p "$RUN_DIR"

    echo "[$(date +%T)] ${RUN_ID} (t=${THREADS} g=${GPUS} trace=${TRACE})"

    EXTRA_FLAGS=()
    case "$runtime" in
        starpu)
            export STARPU_SCHED="$scheduler"
            export STARPU_NCUDA="$GPUS"
            export STARPU_CALIBRATE=0   # freeze the model: measure, don't fit
            unset PARSEC_MCA_mca_sched PARSEC_MCA_device_cuda_enabled
            unset PARSEC_MCA_mca_pins PARSEC_MCA_device_show_statistics CHAMELEON_PARSEC_DOT PARSEC_MCA_profile_filename
            if [[ "$TRACE_STATS" == "1" ]]; then
                export STARPU_BUS_STATS=1 STARPU_WORKER_STATS=1
            else
                unset STARPU_BUS_STATS STARPU_WORKER_STATS
            fi
            if [[ "$TRACE" == "1" ]]; then
                export STARPU_FXT_TRACE=1
                export STARPU_FXT_PREFIX="${RUN_DIR}/"
                EXTRA_FLAGS+=(--trace)
            else
                export STARPU_FXT_TRACE=0
                unset STARPU_FXT_PREFIX
            fi
            ;;
        parsec)
            export PARSEC_MCA_mca_sched="$scheduler"
            export PARSEC_MCA_device_cuda_enabled=1
            unset STARPU_SCHED STARPU_NCUDA STARPU_CALIBRATE
            unset STARPU_BUS_STATS STARPU_WORKER_STATS STARPU_FXT_TRACE STARPU_FXT_PREFIX
            if [[ "$TRACE_STATS" == "1" ]]; then
                export PARSEC_MCA_device_show_statistics=1
            else
                unset PARSEC_MCA_device_show_statistics
            fi
            if [[ "$TRACE" == "1" ]]; then
                export PARSEC_MCA_profile_filename="${RUN_DIR}/cham_${kernel}"
                if [[ "$TRACE_FULL" == "1" ]]; then
                    export PARSEC_MCA_mca_pins=task_profiler
                else
                    unset PARSEC_MCA_mca_pins
                fi
                if [[ "$TRACE_DAG" == "1" ]]; then
                    # Chameleon forwards "-." to parsec_init when this is set, so
                    # the grapher writes ${RUN_DIR}/cham_${kernel}.dot.
                    export CHAMELEON_PARSEC_DOT="${RUN_DIR}/cham_${kernel}"
                else
                    unset CHAMELEON_PARSEC_DOT
                fi
            else
                unset PARSEC_MCA_profile_filename PARSEC_MCA_mca_pins CHAMELEON_PARSEC_DOT
            fi
            ;;
        *)
            echo "unknown runtime: $runtime" >&2
            exit 1
            ;;
    esac

    set +e
    "$bin" -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
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
