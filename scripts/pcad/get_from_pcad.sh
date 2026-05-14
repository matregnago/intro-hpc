#!/bin/bash
#
set -euo pipefail

[[ -f flake.nix ]] || { echo "rode a partir da raiz do projeto"; exit 1; }

SSH_HOST="pcad"
REMOTE_DIR="intro_hpc"

rsync --verbose --progress --recursive --links --times \
    --exclude='.git/' \
    --exclude='.claude/' \
    --exclude='.venv/' \
    --exclude='.direnv/' \
    --exclude='.ipynb_checkpoints/' \
    --exclude='__pycache__/' \
    --exclude='results_local/' \
    "${SSH_HOST}:~/${REMOTE_DIR}/" ./

echo "Dados movidos com sucesso!"
