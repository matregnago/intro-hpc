#!/bin/bash
# Diagnostic: isolate the chameleon_stesting segfault seen in job 796182.
# Small problem so each case finishes in seconds.
set -uo pipefail

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PARSEC_MCA_device_cuda_enabled=1
export PARSEC_MCA_mca_sched=lfq

N=2880; B=480; T=8
BIN=$(command -v chameleon_stesting)
echo "chameleon_stesting = $BIN"
echo "gdb present        = $(command -v gdb || echo NO)"
echo "N=$N B=$B T=$T"
echo

run() {
    local name=$1; shift
    echo "==================== CASE: $name ===================="
    echo "env: CHAMELEON_PARSEC_DOT=${CHAMELEON_PARSEC_DOT-<unset>} pins=${PARSEC_MCA_mca_pins-<unset>} prof=${PARSEC_MCA_profile_filename-<unset>} gpu_arg=$GPUARG"
    stdbuf -oL -eL "$@" -o spotrf -n "$N" -b "$B" -t "$T" -g "$GPUARG" --nowarmup
    echo "---> exit=$? ($name)"
    echo
}

# CASE 1: baseline, GPU on, no DOT  -> expected to SUCCEED (matches prior job 795377)
unset CHAMELEON_PARSEC_DOT PARSEC_MCA_mca_pins PARSEC_MCA_profile_filename
GPUARG=1
run baseline_gpu_noDOT "$BIN"

# CASE 2: DOT + GPU  -> reproduces job 796182
export CHAMELEON_PARSEC_DOT=/tmp/repro_cham_c2
GPUARG=1
run DOT_gpu "$BIN"

# CASE 3: DOT + CPU only -> is it DOT alone, or DOT*GPU interaction?
export CHAMELEON_PARSEC_DOT=/tmp/repro_cham_c3
GPUARG=0
run DOT_cpu "$BIN"

# CASE 4: full capture env (DOT + task_profiler + profile_filename + GPU) -> exact 796182 conditions
export CHAMELEON_PARSEC_DOT=/tmp/repro_cham_c4
export PARSEC_MCA_mca_pins=task_profiler
export PARSEC_MCA_profile_filename=/tmp/repro_prof_c4
GPUARG=1
run DOT_gpu_pins "$BIN"

# CASE 5: backtrace of the failing case, if gdb available
if command -v gdb >/dev/null; then
    echo "==================== CASE: gdb backtrace (DOT+GPU) ===================="
    export CHAMELEON_PARSEC_DOT=/tmp/repro_cham_c5
    unset PARSEC_MCA_mca_pins PARSEC_MCA_profile_filename
    gdb -batch -nx \
        -ex 'set pagination off' \
        -ex run \
        -ex 'bt' \
        -ex 'info registers' \
        --args "$BIN" -o spotrf -n "$N" -b "$B" -t "$T" -g 1 --nowarmup
    echo "---> gdb done"
fi
echo "ALL CASES DONE"
