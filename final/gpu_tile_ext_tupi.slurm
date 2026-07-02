#!/bin/bash -l
#SBATCH --partition=tupi
#SBATCH --time=08:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --job-name=gpu_tile_ext_tupi
#SBATCH --nodelist=tupi3,tupi4,tupi5,tupi6

# E1b -- varredura de tile ESTENDIDA (N=40000) na GPU, particao TUPI
# (i9-14900KF 24c, RTX 4090 24GB). THREADS=23 (reserva 1 core p/ dirigir a GPU).
#
# Roda final/doe_gpu_tile_ext_tupi.csv (480 runs: 4 ops x 4 sched x 10 tiles x
# 3 reps) nos dois runtimes, SEM trace -- MESMO grid/design da particao poti
# (final/gpu_tile_ext_poti.slurm), so muda o no e THREADS, p/ a comparacao
# 4070 vs 4090 ser no mesmo problema. Continua o E1 (N=20000): aqui o N maior
# mantem a grade de tiles saudavel ate b=4000 (10x10), isolando o efeito do
# tile do efeito de perda de paralelismo. Os modelos de desempenho do StarPU
# (dmda/dmdas) sao CALIBRADOS por (precisao,algoritmo,n,b) antes das medidas
# (final/gpu_doe_sweep.sh::calibrate_starpu, modelo isolado em
# $RESULTS_DIR/.starpu, congelado p/ as runs cronometradas).
#
# N=40000 e ~8x mais caro por run que o E1 (N^3) e dobra o nro de runs -> tempo
# de parede maior; por isso --time=08:00:00 (vs 04:00:00 do E1).
#
# OBS: o run anterior da tupi (797982) teve a GPU aparentemente disputada
# (StarPU OOM na calibracao de banda + GFlops do PaRSEC abaixo do 4070); confira
# que a GPU do no esta livre antes de submeter.
#
#   sbatch final/gpu_tile_ext_tupi.slurm

set -euo pipefail

NODE=tupi
THREADS=23
GPUS=1
DESIGN_SRC="final/doe_gpu_tile_ext_${NODE}.csv"

echo "==============================================================="
echo "Date        : $(date)"
echo "JobID       : $SLURM_JOB_ID   Partition: $SLURM_JOB_PARTITION"
echo "Node(s)     : $SLURM_JOB_NODELIST   THREADS=$THREADS GPUS=$GPUS"
echo "==============================================================="

WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_RESULTS_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

[[ -f "$DESIGN_SRC" ]] || { echo "DoE nao encontrado: $DESIGN_SRC" >&2; exit 1; }

copyback() {
    echo "$LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    [ -d "$LOCAL_RESULTS_DIR" ] && cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
}
trap copyback EXIT

export RESULTS_DIR="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
cp "$DESIGN_SRC" "${RESULTS_DIR}/design.csv"
export DESIGN_FILE="${RESULTS_DIR}/design.csv"

# E1b nao coleta trace -- so time/gflops.
export TRACE=0 TRACE_FULL=0 TRACE_DAG=0 TRACE_STATS=0
export THREADS GPUS

run_phase() {
    local shell="$1" runtime="$2"
    echo ">>> fase: ${runtime}  (nixw nix develop .#${shell})"
    if ! nixw nix develop ".#${shell}" --impure --command env \
        RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
        THREADS="$THREADS" GPUS="$GPUS" \
        TRACE="$TRACE" TRACE_FULL="$TRACE_FULL" TRACE_DAG="$TRACE_DAG" \
        TRACE_STATS="$TRACE_STATS" \
        bash final/gpu_doe_sweep.sh "$runtime"; then
        echo "!!! fase ${runtime} (.#${shell}) falhou; continuando" >&2
    fi
}

# StarPU primeiro (calibra dmda/dmdas), depois PaRSEC.
run_phase starpu starpu
run_phase parsec parsec

echo "==> E1b (tupi) concluido."
