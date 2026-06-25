#!/bin/bash -l
#SBATCH --partition=poti
#SBATCH --time=16:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --job-name=chameleon-gpu-phases
# poti ja expoe o RTX 4070 ao job; se a particao exigir alocacao explicita de
# GPU, descomente a linha abaixo:
##SBATCH --gres=gpu:1

set -euo pipefail

echo "==============================================================="
echo "Date        : $(date)"
echo "SubmitHost  : $(hostname)"
echo "SubmitDir   : $SLURM_SUBMIT_DIR"
echo "JobID       : $SLURM_JOB_ID"
echo "Node(s)     : $SLURM_JOB_NODELIST"
echo "==============================================================="

# Roda as etapas E1-E4 do plano de experimentos GPU (docs/plano_experimentos_gpu.md)
# no mesmo molde dos jobs existentes (slurm/traces.slurm / cuda_vs_cpu.slurm):
# rsync do submit dir -> $SCRATCH, entra nos dev shells CUDA via `nixw nix
# develop`, e um trap EXIT copia os resultados de volta. A E4 roda com TRACE=1
# para capturar os rastros (StarPU FxT / PaRSEC profile) -- a peca que faltava
# no gpu_doe_sweep.sh e que fazia a captura local nao funcionar.

WORK_DIR="$SCRATCH/intro_hpc"
LOCAL_RESULTS_DIR="$WORK_DIR/results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
FINAL_RESULTS_DIR="$SLURM_SUBMIT_DIR/data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$LOCAL_RESULTS_DIR" "$FINAL_RESULTS_DIR"
cp -r "$SLURM_SUBMIT_DIR/." "$WORK_DIR/"
cd "$WORK_DIR"

copyback() {
    echo "$LOCAL_RESULTS_DIR -> $FINAL_RESULTS_DIR"
    if [ -d "$LOCAL_RESULTS_DIR" ]; then
        cp -r "$LOCAL_RESULTS_DIR/." "$FINAL_RESULTS_DIR/" || true
    fi
}
trap copyback EXIT

# Knobs compartilhados (override aqui). 16 threads CPU + 1 GPU espelham os runs
# documentados no poti (docs/parsec_cuda_results.md).
export THREADS="${THREADS:-16}"
export GPUS="${GPUS:-1}"

# Escalabilidade em N (E3): FP32 (onde a GPU de consumo ganha da CPU). Ns sao
# multiplos de 960 -> tiling limpo com b=960.
E3_BIN="${E3_BIN:-chameleon_stesting}"
E3_OPS="${E3_OPS:-spotrf sgeqrf}"
E3_NS="${E3_NS:-9600 19200 28800 38400 48000}"
E3_B="${E3_B:-960}"
E3_REPS="${E3_REPS:-3}"

RES_ROOT="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

# ---------------------------------------------------------------------------
# Fase DoE-driven (E1/E2/E4): roda gpu_doe_sweep.sh para cada runtime dentro do
# seu dev shell CUDA. As fases sao independentes -- uma falha (ex.: build
# faltando) e logada mas nao aborta o resto, entao uma submissao ainda gera o
# que deu certo. Cada nix develop e um processo filho, entao STARPU_HOME (o
# cache de modelos do dmda/dmdas) e recalculado fresco por fase, sem poluir.
#   $1 = nome da fase   $2 = design csv (relativo ao WORK_DIR)   $3 = TRACE (0/1)
#   $4 = TRACE_FULL (0/1, opcional): adiciona PINS task_profiler + DAG (.dot) p/ StarVZ
run_doe_phase() {
    local name="$1" design_src="$2" trace="$3" trace_full="${4:-0}"
    local rdir="${RES_ROOT}/${name}"

    if [[ ! -f "$design_src" ]]; then
        echo "!!! pulando ${name}: design '${design_src}' nao encontrado" >&2
        echo "    gere com o Rscript correspondente antes de sincronizar (ver docs)" >&2
        return 0
    fi

    mkdir -p "$rdir/runs"
    cp "$design_src" "$rdir/design.csv"

    export RESULTS_DIR="$rdir"
    export RESULTS_FILE="$rdir/results.csv"
    export RUNS_DIR="$rdir/runs"
    export DESIGN_FILE="$rdir/design.csv"
    export TRACE="$trace"
    export TRACE_FULL="$trace_full"

    echo
    echo "########## fase ${name} (trace=${trace} trace_full=${trace_full}) ##########"
    for combo in "starpu-cuda:starpu" "parsec:parsec"; do
        local shell="${combo%%:*}" rt="${combo##*:}"
        echo ">>> ${name}: nix develop .#${shell} -> gpu_doe_sweep.sh ${rt}"
        if ! nixw nix develop ".#${shell}" --impure --command bash -c "
            set -euo pipefail
            bash scripts/gpu_doe_sweep.sh ${rt}
        "; then
            echo "!!! fase ${name}/${rt} (.#${shell}) falhou; continuando" >&2
        fi
    done
}

# ---------------------------------------------------------------------------
# Fase de escalabilidade (E3): gpu_sweep.sh (sem DoE CSV, varre N por env).
run_e3_phase() {
    local name="e3_scal"
    local rdir="${RES_ROOT}/${name}"
    mkdir -p "$rdir/runs"

    export RESULTS_DIR="$rdir"
    export RESULTS_FILE="$rdir/results.csv"
    export RUNS_DIR="$rdir/runs"
    export TESTING_BIN="$E3_BIN" OPS="$E3_OPS" NS="$E3_NS" B="$E3_B" REPS="$E3_REPS"
    unset DESIGN_FILE TRACE || true

    echo
    echo "########## fase ${name} (escalabilidade em N) ##########"
    for combo in "starpu-cuda:starpu" "parsec:parsec"; do
        local shell="${combo%%:*}" rt="${combo##*:}"
        echo ">>> ${name}: nix develop .#${shell} -> gpu_sweep.sh ${rt}"
        if ! nixw nix develop ".#${shell}" --impure --command bash -c "
            set -euo pipefail
            bash scripts/gpu_sweep.sh ${rt}
        "; then
            echo "!!! fase ${name}/${rt} (.#${shell}) falhou; continuando" >&2
        fi
    done
}

# E1 -- selecao de tile (b)
run_doe_phase  e1_tile    "scripts/doe/doe_gpu_tile.csv"    0
# E2 -- varredura de escalonadores (figura-mae)
run_doe_phase  e2_sched   "scripts/doe/doe_gpu_sched.csv"   0
# E3 -- escalabilidade em N
run_e3_phase
# E4 -- rastros (StarVZ) com TRACE=1 e TRACE_FULL=1 (PINS + DAG .dot do PaRSEC)
run_doe_phase  e4_traces  "scripts/doe/doe_gpu_traces.csv"  1 1

echo
echo "==> Fim. Resultados por fase em data/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/{e1_tile,e2_sched,e3_scal,e4_traces}/"
