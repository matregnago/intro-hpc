#!/bin/bash

# Roda NO FRONT-END do PCAD (gppd-hpc). Descobre o no da sua alocacao
# interativa (o job 'sleep infinity' na tupi) e dispara nele
# scripts/pcad/collect_n_size_tupi_node.sh, que roda o experimento E2b
# (escalonamento vs tamanho do problema N, no pico de tile b=1000 de E1 tupi) no
# $SCRATCH e copia de volta pra ~/intro_hpc/data/.
#
#   uso: bash scripts/pcad/collect_n_size_tupi_pcad.sh <JOBID>
#
# JOBID e obrigatorio (passe o do 'sleep infinity' na tupi).
# DOE_CSV (opcional): CSV de design a rodar (default final/doe_gpu_n_size_tupi.csv).
#   Use para rodar so um subset, p.ex. complementar Ns de um job anterior:
#   DOE_CSV=final/doe_gpu_n_size_tupi_extend.csv bash ...pcad.sh <JOBID>
#
# Depois, do SEU PC LOCAL: bash scripts/pcad/get_from_pcad.sh

set -euo pipefail

SUBMIT_DIR="${SUBMIT_DIR:-$HOME/intro_hpc}"
MECH="${MECH:-srun}"
JOBTAG="${JOBTAG:-$(date +%Y%m%d_%H%M%S)}"
GPUS="${GPUS:-1}"
DOE_CSV="${DOE_CSV:-final/doe_gpu_n_size_tupi.csv}"

NODE_SCRIPT="scripts/pcad/collect_n_size_tupi_node.sh"
[[ -f "$SUBMIT_DIR/$NODE_SCRIPT" ]] || {
    echo "nao achei $SUBMIT_DIR/$NODE_SCRIPT" >&2
    echo "sincronize o repo: bash scripts/pcad/sync_to_pcad.sh (do seu PC local)" >&2
    exit 1
}

# Descobre o no a partir do JOBID.
JOBID="${1:?uso: $0 <JOBID>}"
NODE="$(squeue -h -j "$JOBID" -o '%N')"
[[ -n "$NODE" ]] || { echo "job $JOBID nao encontrado (squeue vazio)" >&2; exit 1; }

echo "==> Job $JOBID  no(s): ${NODE:-?}  (mech=$MECH, tag=$JOBTAG)"

case "$MECH" in
    srun)
        # --overlap: compartilha os recursos com o passo 'sleep infinity' ativo.
        srun --jobid "$JOBID" --overlap --chdir="$SUBMIT_DIR" \
            env JOBTAG="$JOBTAG" SUBMIT_DIR="$SUBMIT_DIR" \
                GPUS="$GPUS" DOE_CSV="$DOE_CSV" \
            bash "$NODE_SCRIPT"
        ;;
    ssh)
        node1="${NODE%%,*}"   # primeiro no, caso seja uma lista
        [[ -n "$node1" ]] || { echo "nao consegui resolver o no via squeue" >&2; exit 1; }
        ssh "$node1" \
            "JOBTAG='$JOBTAG' SUBMIT_DIR='$SUBMIT_DIR' \
             SCRATCH=\"\${SCRATCH:-\$HOME/scratch}\" \
             GPUS='$GPUS' DOE_CSV='$DOE_CSV' \
             bash '$SUBMIT_DIR/$NODE_SCRIPT'"
        ;;
    *) echo "MECH invalido: $MECH (use srun|ssh)" >&2; exit 1 ;;
esac

echo
echo "==> Pronto. Resultados em $SUBMIT_DIR/data/gpu_n_size_tupi_$JOBTAG"
echo "    Do seu PC local, puxe com: bash scripts/pcad/get_from_pcad.sh"
