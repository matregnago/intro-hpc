#!/usr/bin/env bash
# Pre-populate the local Nix store with tarballs from the shared filesystem.
# This avoids network downloads on compute nodes that have no internet access.
# Called from SLURM scripts before any `nix develop` or `nix build` commands.
#
# Usage: nix-prefetch-tarballs.sh <tarballs-dir>
set -euo pipefail

TARBALLS_DIR="${1:?Usage: $0 <tarballs-dir>}"

if [ ! -d "$TARBALLS_DIR" ]; then
  echo "Warning: tarballs directory not found: $TARBALLS_DIR" >&2
  echo "Run scripts/download-tarballs.sh on the login node first." >&2
  exit 1
fi

for f in "$TARBALLS_DIR"/*.tar.gz "$TARBALLS_DIR"/*.tar.bz2 "$TARBALLS_DIR"/*.tar.xz; do
  [ -f "$f" ] || continue
  name="$(basename "$f")"
  echo "Adding $name to Nix store ..."
  nixw nix store add --hash-algo sha256 --mode flat --name "$name" "$f"
done

echo "Nix store pre-populated from $TARBALLS_DIR"
