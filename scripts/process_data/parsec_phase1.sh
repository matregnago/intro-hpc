#!/usr/bin/env bash
#
# Offline phase 1 for PaRSEC runs: reconstruct the StarVZ parquet set directly
# from cham_*.dot + cham_*.prof-* via dbp2xml. Use when the official phase 1
# can't produce the tables, i.e. the GPU runs (dbp2paje breaks on the GPU paje).
#
# Per run_dir, in order:
#   parsec_tasks_to_parquet.r        -> tasks.parquet + states.parquet
#   parsec_dag_to_parquet.r          -> dag.parquet
#   parsec_application_to_parquet.r  -> application.parquet + y/colors.parquet
#   parsec_ready_to_parquet.r        -> variable.parquet (reconstructed Ready)
#
# dbp2xml is resolved once and exported, so the converters reuse it.
#
# Usage:
#   scripts/process_data/parsec_phase1.sh <run_dir> [run_dir...]

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <run_dir> [run_dir...]" >&2
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

# Resolve dbp2xml: $DBP2XML, PATH, or build the flake's .#parsec.
# --impure: parsec.nix reads its source from $PWD.
if [[ -z "${DBP2XML:-}" || ! -x "${DBP2XML:-}" ]]; then
    if command -v dbp2xml >/dev/null 2>&1; then
        DBP2XML="$(command -v dbp2xml)"
    else
        echo "dbp2xml not on PATH; building .#parsec (first run only)..." >&2
        out="$(cd "$REPO" && PWD="$REPO" nix build --impure --no-link \
            --print-out-paths ".#parsec" 2>/dev/null | tail -1)"
        DBP2XML="${out:+$out/bin/dbp2xml}"
    fi
fi
if [[ -z "${DBP2XML:-}" || ! -x "$DBP2XML" ]]; then
    echo "ERROR: could not obtain dbp2xml (put it on PATH, set \$DBP2XML, or run" \
         "from a tree where 'nix build --impure .#parsec' works)" >&2
    exit 1
fi
export DBP2XML
echo "dbp2xml: $DBP2XML"

phase1_one() {
    local run_dir="${1%/}"
    if [[ ! -d "$run_dir" ]]; then
        echo "!!! skipping (not a dir): $run_dir" >&2
        return 1
    fi
    shopt -s nullglob
    local dot=("$run_dir"/cham_*.dot)
    local prof=("$run_dir"/cham_*.prof-*)
    shopt -u nullglob
    # the converters expect an unpacked .prof; ignore a compressed-only one
    local p clean=()
    for p in "${prof[@]}"; do [[ "$p" == *.tar.gz ]] || clean+=("$p"); done
    prof=("${clean[@]}")

    if [[ ${#dot[@]} -eq 0 || ${#prof[@]} -eq 0 ]]; then
        echo "!!! skipping (need cham_*.dot AND cham_*.prof-*): $run_dir" >&2
        return 1
    fi

    echo "==> $run_dir"
    Rscript "$HERE/parsec_tasks_to_parquet.r"       "$run_dir"
    Rscript "$HERE/parsec_dag_to_parquet.r"         "$run_dir"
    Rscript "$HERE/parsec_application_to_parquet.r" "$run_dir"
    Rscript "$HERE/parsec_ready_to_parquet.r"       "$run_dir"
    echo "    done: $(ls "$run_dir"/*.parquet 2>/dev/null | wc -l) parquet files"
}

rc=0
for d in "$@"; do
    phase1_one "$d" || rc=1
done
exit "$rc"
