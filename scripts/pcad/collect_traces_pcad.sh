#!/bin/bash

# Roda NO FRONT-END do PCAD (gppd-hpc). Descobre o no da sua alocacao
# interativa (o job 'sleep infinity') e dispara nele
# scripts/pcad/collect_traces_node.sh, que coleta os traces no $SCRATCH e
# copia de volta pra ~/intro_hpc/data/ (o $HOME e compartilhado, entao os
# resultados ficam visiveis aqui no front-end).
#
#   uso: bash scripts/pcad/collect_traces_pcad.sh [JOBID]
#
# Sem JOBID, pega seu unico job RUNNING. Mecanismo via env MECH:
#   srun (default) -- srun --jobid --overlap: roda como um passo dentro da
#                     alocacao, herdando \$SCRATCH e o binding da GPU. Mais robusto.
#   ssh            -- ssh direto no no (so se srun --overlap nao funcionar aqui).
#
# Depois, do SEU PC LOCAL: bash scripts/pcad/get_from_pcad.sh

set -euo pipefail

SUBMIT_DIR="${SUBMIT_DIR:-$HOME/intro_hpc}"
MECH="${MECH:-srun}"
JOBTAG="${JOBTAG:-$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"
PHASE="${PHASE:-both}"

NODE_SCRIPT="scripts/pcad/collect_traces_node.sh"
[[ -f "$SUBMIT_DIR/$NODE_SCRIPT" ]] || {
    echo "nao achei $SUBMIT_DIR/$NODE_SCRIPT" >&2
    echo "sincronize o repo: bash scripts/pcad/sync_to_pcad.sh (do seu PC local)" >&2
    exit 1
}

# Descobre JOBID + no(s).
JOBID="${1:-}"
if [[ -z "$JOBID" ]]; then
    mapfile -t rows < <(squeue --me -h -t RUNNING -o '%i|%j|%N')
    if [[ "${#rows[@]}" -eq 0 ]]; then
        echo "nenhum job RUNNING seu. Lance antes: sbatch -p poti --wrap 'sleep infinity'" >&2
        exit 1
    elif [[ "${#rows[@]}" -gt 1 ]]; then
        echo "varios jobs RUNNING; passe o JOBID desejado:" >&2
        printf '  %s\n' "${rows[@]}" >&2
        exit 1
    fi
    IFS='|' read -r JOBID JOBNAME NODE <<< "${rows[0]}"
else
    NODE="$(squeue -h -j "$JOBID" -o '%N')"
fi

echo "==> Job $JOBID  no(s): ${NODE:-?}  (mech=$MECH, tag=$JOBTAG)"
echo "==> Sairao em: $SUBMIT_DIR/data/manual_traces_$JOBTAG"

case "$MECH" in
    srun)
        # --overlap: compartilha os recursos com o passo 'sleep infinity' ativo.
        srun --jobid "$JOBID" --overlap --chdir="$SUBMIT_DIR" \
            env JOBTAG="$JOBTAG" SUBMIT_DIR="$SUBMIT_DIR" \
                THREADS="$THREADS" GPUS="$GPUS" PHASE="$PHASE" \
            bash "$NODE_SCRIPT"
        ;;
    ssh)
        node1="${NODE%%,*}"   # primeiro no, caso seja uma lista
        [[ -n "$node1" ]] || { echo "nao consegui resolver o no via squeue" >&2; exit 1; }
        # SCRATCH e avaliado NO NO (se o profile dele definir, usa; senao fallback).
        ssh "$node1" \
            "JOBTAG='$JOBTAG' SUBMIT_DIR='$SUBMIT_DIR' \
             SCRATCH=\"\${SCRATCH:-\$HOME/scratch}\" \
             THREADS='$THREADS' GPUS='$GPUS' PHASE='$PHASE' \
             bash '$SUBMIT_DIR/$NODE_SCRIPT'"
        ;;
    *) echo "MECH invalido: $MECH (use srun|ssh)" >&2; exit 1 ;;
esac

echo
echo "==> Pronto. Traces em $SUBMIT_DIR/data/manual_traces_$JOBTAG"
echo "    Do seu PC local, puxe com: bash scripts/pcad/get_from_pcad.sh"
