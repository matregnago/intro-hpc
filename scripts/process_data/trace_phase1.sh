#!/usr/bin/env bash
#
# StarVZ phase 1 for Chameleon trace runs (StarPU or PaRSEC, auto-detected):
# writes the phase-1 *.parquet files next to each run's raw trace.
#   PaRSEC -> dbp2paje converts cham_*.prof-* to paje.trace, then `starvz -1 -t`.
#   StarPU -> `starvz -1` runs starpu_fxt_tool on the prof_file_* itself.
#
# dbp2paje MUST get -n (--name-all-containers): without it the worker/stream
# containers get empty aliases, every state ends up with ResourceId " ", and
# phase 1 dies in hl_y_paje_tree.
#
# Usage:
#   scripts/process_data/trace_phase1.sh <run_dir> [run_dir...]

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <run_dir> [run_dir...]" >&2
    exit 1
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Resolve a tool: prefer PATH, otherwise build the flake package and point at
# its binary. --impure: parsec.nix reads its source from $PWD.
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
