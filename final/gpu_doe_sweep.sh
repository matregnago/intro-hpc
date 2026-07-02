#!/bin/bash

# DoE-driven GPU sweep for the StarPU-vs-PaRSEC scheduler study. Run *inside*
# the matching dev shell (.#starpu or .#parsec). Reads $DESIGN_FILE (a CSV
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
# Tracing for E4. TRACE=1 captures StarPU FxT / PaRSEC profile files next to
# each run's log -- exactly the mechanism scripts/run.sh uses -- so that
# scripts/analysis/plot_traces.r (StarVZ) can read them. Default off.
TRACE="${TRACE:-0}"
# TRACE_FULL=1 (only meaningful with TRACE=1) collects the *richest* trace we can
# without recompiling, for the StarVZ "Fase A" adapter (docs/starvz_e_captura_parsec.md).
# For PaRSEC it adds the task_profiler PINS module, so the .prof also carries
# runtime/scheduler states (EXEC/RELEASE_DEPS/ACTIVATE_CB) -> a non-empty
# `starpu`-equivalent tibble. CUDA H2D/D2H/exec/own are already all-on by default
# (parsec/devices/cuda/dev_cuda.c: parsec_cuda_trackable_events). StarPU's FxT is
# already complete under TRACE=1, so no extra knob is needed on that side.
# It also collects the PaRSEC task-graph DAG (.dot) -- the old "Fase B": the
# build already enables PARSEC_PROF_GRAPHER (parsec.nix) and Chameleon forwards
# the "-." option to parsec_init when CHAMELEON_PARSEC_DOT is set
# (runtime/parsec/control/runtime_control.c), so the grapher writes
# <run_dir>/cham_<kernel>.dot next to the .prof -- the piece StarVZ needs for
# the Dag/Task_handles tibbles. Default off.
TRACE_FULL="${TRACE_FULL:-0}"
# TRACE_DAG=1 (only meaningful with TRACE=1, PaRSEC only) emits the task-graph
# DAG (.dot) via CHAMELEON_PARSEC_DOT. Decoupled from TRACE_FULL because at
# experiment scale (n=19200) the DAG is a huge unrenderable hairball and the
# grapher pollutes the timed trace -- the DAG is meant for a small dedicated run
# (scripts/capture_parsec_dag.sh / slurm/parsec_dag.slurm). Defaults to
# TRACE_FULL so existing callers keep their behavior.
TRACE_DAG="${TRACE_DAG:-$TRACE_FULL}"
# TRACE_STATS=1 prints the runtimes' NATIVE data-movement counters into each run's
# log (the cleanest H2D/D2H evidence for the GPU thrashing analysis, independent
# of the trace): StarPU STARPU_BUS_STATS (+ WORKER_STATS), PaRSEC
# device_show_statistics. They dump at deinit to stderr/stdout (captured in the
# log) and add no measurable overhead, so default them ON whenever tracing.
TRACE_STATS="${TRACE_STATS:-$TRACE}"
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

# Validated env (local-gpu-compare memory): one BLAS thread per worker.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export STARPU_SILENT=1

# Workers (CUDA streams) per GPU. StarPU defaults to 1 execution stream per
# device; PaRSEC defaults to 4 (device_cuda_max_streams=6 -> 2 for transfer + 4
# for kernels), so match it at 4 for a paired CPU+GPU comparison -- otherwise
# StarPU shows a single GPU lane vs PaRSEC's 4 and the StarVZ space-time panels
# aren't comparable. Exported here (not in the per-run branch) so the dmda/dmdas
# CALIBRATION sees the same worker config as the timed runs. No-op on CPU-only
# runs (no CUDA device); ignored by PaRSEC. Override via env.
export STARPU_NWORKER_PER_CUDA="${STARPU_NWORKER_PER_CUDA:-4}"

# CPU worker binding. On native Linux (PCAD poti/tupi) we WANT StarPU to bind
# workers to cores for deterministic placement -- so binding stays ON by
# default. Under WSL2 hwloc misdetects the topology, which makes binding crash/
# mislay workers (see parsec-cuda-trace-local), so we disable it THERE only.
# Auto-detect WSL; an explicit STARPU_WORKERS_NOBIND in the environment wins.
if [[ -z "${STARPU_WORKERS_NOBIND:-}" ]] && grep -qiE 'microsoft|wsl' /proc/version 2>/dev/null; then
    export STARPU_WORKERS_NOBIND=1
    export STARPU_WORKERS_GETBIND=0
fi

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

    EXTRA_FLAGS=()
    case "$runtime" in
        starpu)
            export STARPU_SCHED="$scheduler"
            export STARPU_NCUDA="$GPUS"
            export STARPU_CALIBRATE=0   # freeze the model: measure, don't fit
            unset PARSEC_MCA_mca_sched PARSEC_MCA_device_cuda_enabled || true
            if [[ "$TRACE_STATS" == "1" ]]; then
                export STARPU_BUS_STATS=1 STARPU_WORKER_STATS=1
            else
                unset STARPU_BUS_STATS STARPU_WORKER_STATS || true
            fi
            if [[ "$TRACE" == "1" ]]; then
                export STARPU_FXT_TRACE=1
                export STARPU_FXT_PREFIX="${RUN_DIR}/"
                EXTRA_FLAGS+=(--trace)
            else
                export STARPU_FXT_TRACE=0
                unset STARPU_FXT_PREFIX || true
            fi
            ;;
        parsec)
            export PARSEC_MCA_mca_sched="$scheduler"
            export PARSEC_MCA_device_cuda_enabled=1
            unset STARPU_SCHED STARPU_NCUDA STARPU_CALIBRATE || true
            if [[ "$TRACE_STATS" == "1" ]]; then
                export PARSEC_MCA_device_show_statistics=1
            else
                unset PARSEC_MCA_device_show_statistics || true
            fi
            if [[ "$TRACE" == "1" ]]; then
                export PARSEC_MCA_profile_filename="${RUN_DIR}/cham_${kernel}"
                if [[ "$TRACE_FULL" == "1" ]]; then
                    export PARSEC_MCA_mca_pins=task_profiler
                else
                    unset PARSEC_MCA_mca_pins || true
                fi
                if [[ "$TRACE_DAG" == "1" ]]; then
                    # DAG (.dot): Chameleon forwards "-." to parsec_init when this
                    # is set, so the grapher writes ${RUN_DIR}/cham_${kernel}.dot.
                    export CHAMELEON_PARSEC_DOT="${RUN_DIR}/cham_${kernel}"
                else
                    unset CHAMELEON_PARSEC_DOT || true
                fi
            else
                unset PARSEC_MCA_profile_filename PARSEC_MCA_mca_pins CHAMELEON_PARSEC_DOT || true
            fi
            ;;
        *) echo "unknown runtime: $runtime" >&2; exit 1 ;;
    esac

    echo "[$(date +%T)] ${RUN_ID} (t=${THREADS} g=${GPUS} trace=${TRACE})"
    set +e
    "$bin" -o "$kernel" -n "$n" -b "$b" -t "$THREADS" -g "$GPUS" --nowarmup \
        "${EXTRA_FLAGS[@]}" > "$LOG_FILE" 2>&1
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
