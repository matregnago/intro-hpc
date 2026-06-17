#!/usr/bin/env bash
#
# StarVZ phase 1 for Chameleon trace runs (PaRSEC or StarPU, auto-detected).
#
# Each run dir holds a raw trace captured by chameleon_dtesting with TRACE=1
# (scripts/gpu_doe_sweep.sh):
#   - PaRSEC: a binary profile  cham_<kernel>-0.prof-XXXXXX
#   - StarPU: one or more FxT files  prof_file_*
#
# StarVZ cannot read either format directly. This script runs phase 1 and writes
# the phase-1 *.parquet files next to the trace (application.parquet,
# starpu.parquet, y.parquet, ...), which starvz_read() and
# scripts/analysis/plot_trace.r consume.
#
#   PaRSEC  -> dbp2paje converts the .prof to paje.trace, then
#              `starvz -1 -t` (--use-paje-trace) runs phase 1 on it.
#   StarPU  -> `starvz -1` runs starpu_fxt_tool on the prof_file_* itself.
#
# IMPORTANT (PaRSEC): dbp2paje is invoked with -n (--name-all-containers).
# Without it the worker/stream containers get empty aliases, parsec_paje_fix.sh
# dedups them into one, every state ends up with ResourceId " ", and phase 1
# dies in hl_y_paje_tree ("Column `Nature` doesn't exist"). -n keeps them apart.
#
# Usage:
#   scripts/analysis/trace_phase1.sh <run_dir> [run_dir...]
#   scripts/analysis/trace_phase1.sh data/chameleon-gpu-phases_795377/e4_traces/runs/00*_*dpotrf*
#
# dbp2paje comes from the `parsec-cpu` flake package and `starvz` from
# `starvz-tools`; both are resolved automatically (built on first use, then
# cached) if not already on PATH, so this also works from a plain `nix develop`.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <run_dir> [run_dir...]" >&2
    exit 1
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Resolve a tool: prefer PATH, otherwise build the flake package and point at its
# binary. --impure is required because parsec.nix reads its source from $PWD.
resolve_tool() {
    local bin="$1" attr="$2" out
    if command -v "$bin" >/dev/null 2>&1; then
        command -v "$bin"
        return
    fi
    out="$(cd "$REPO" && PWD="$REPO" nix build --impure --no-link \
        --print-out-paths ".#${attr}" 2>/dev/null | tail -1)"
    if [[ -z "$out" || ! -x "$out/bin/$bin" ]]; then
        echo "ERROR: could not obtain '$bin' (flake package .#${attr})" >&2
        exit 1
    fi
    echo "$out/bin/$bin"
}

STARVZ="$(resolve_tool starvz starvz-tools)"
DBP2PAJE=""  # resolved lazily, only if a PaRSEC run is found
echo "starvz: $STARVZ"

# PaRSEC: cham_<kernel>-0.prof-XXXX -> paje.trace -> starvz -1 -t
phase1_parsec() {
    local run_dir="$1"
    shopt -s nullglob
    local profs=("$run_dir"/cham_*.prof-*)
    if [[ ${#profs[@]} -eq 0 ]]; then
        local tgz=("$run_dir"/cham_*.prof-*.tar.gz)
        if [[ ${#tgz[@]} -gt 0 ]]; then
            echo "    extracting $(basename "${tgz[0]}")"
            tar -xzf "${tgz[0]}" -C "$run_dir"
            profs=("$run_dir"/cham_*.prof-*)
        fi
    fi
    local p clean=()
    for p in "${profs[@]}"; do [[ "$p" == *.tar.gz ]] || clean+=("$p"); done
    profs=("${clean[@]}")
    shopt -u nullglob

    if [[ -z "$DBP2PAJE" ]]; then DBP2PAJE="$(resolve_tool dbp2paje parsec-cpu)"; fi
    echo "    runtime: PaRSEC ; prof: ${profs[*]##*/}"

    local prof_names=()
    for p in "${profs[@]}"; do prof_names+=("$(basename "$p")"); done
    ( cd "$run_dir" && rm -f paje.trace && \
        "$DBP2PAJE" -n -p -s -o paje "${prof_names[@]}" >/dev/null )
    "$STARVZ" -1 -t "$run_dir"
}

# StarPU: prof_file_* (FxT) -> starvz -1 runs starpu_fxt_tool itself
phase1_starpu() {
    local run_dir="$1"
    echo "    runtime: StarPU ; fxt: $(cd "$run_dir" && echo prof_file_*)"
    "$STARVZ" -1 "$run_dir"
}

phase1_one() {
    local run_dir="${1%/}"
    if [[ ! -d "$run_dir" ]]; then
        echo "!!! skipping (not a dir): $run_dir" >&2
        return 1
    fi
    echo "==> $run_dir"

    shopt -s nullglob
    local has_parsec=("$run_dir"/cham_*.prof-*)
    local has_starpu=("$run_dir"/prof_file_*)
    shopt -u nullglob

    if [[ ${#has_parsec[@]} -gt 0 ]]; then
        phase1_parsec "$run_dir"
    elif [[ ${#has_starpu[@]} -gt 0 ]]; then
        phase1_starpu "$run_dir"
    else
        echo "!!! no cham_*.prof-* nor prof_file_* in $run_dir, skipping" >&2
        return 1
    fi

    echo "    done: $(ls "$run_dir"/*.parquet 2>/dev/null | wc -l) parquet files"
}

rc=0
for d in "$@"; do
    phase1_one "$d" || rc=1
done
exit "$rc"
