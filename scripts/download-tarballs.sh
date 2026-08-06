#!/usr/bin/env bash
# Download large tarballs that cannot be committed to git.
# Run this once on the login node (which has internet access) before submitting jobs.
# Tarballs are stored in ./tarballs/ (gitignored) on the shared filesystem.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARBALLS_DIR="${1:-$SCRIPT_DIR/../tarballs}"
mkdir -p "$TARBALLS_DIR"

declare -A TARBALLS=(
  ["starpu-1.4.7.tar.gz"]="http://files.inria.fr/starpu/starpu-1.4.7/starpu-1.4.7.tar.gz"
)

for name in "${!TARBALLS[@]}"; do
  dest="$TARBALLS_DIR/$name"
  if [ -f "$dest" ]; then
    echo "Already present: $dest"
  else
    echo "Downloading $name ..."
    curl -L --fail -o "$dest" "${TARBALLS[$name]}"
    echo "Downloaded: $dest"
  fi
done

echo "All tarballs ready in $TARBALLS_DIR"
