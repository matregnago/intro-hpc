#!/bin/bash

set -euo pipefail

[[ -f flake.nix ]] || { echo "rode a partir da raiz do projeto"; exit 1; }

SSH_HOST="pcad"
REMOTE_DIR="intro_hpc"

rsync --verbose --progress --recursive --links --times \
    --exclude='.git/' \
    --exclude='.claude/' \
    --exclude='.venv/' \
    --exclude='.direnv/' \
    --exclude='traces_parsec' \
    --exclude='full_factorial_780418/' \
    --exclude='results_pcad/' \
    --exclude='tex/' \
    --exclude='slides/' \
    --exclude='analysis.ipynb' \
    --exclude='*.out' \
    --exclude='*.err' \
    ./ "${SSH_HOST}:~/${REMOTE_DIR}/"

echo "Dados movidos com sucesso!"
