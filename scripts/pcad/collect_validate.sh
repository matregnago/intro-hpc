#!/bin/bash
#
# Coleta MINIMA de validacao -- roda interativamente num no da PCAD (poti/RTX 4070)
# para capturar, num caso PEQUENO e rapido, TODOS os dados que as visualizacoes
# precisam, dos dois runtimes. Quando a gente confirmar que o pipeline inteiro
# roda em cima desses traces, voce dispara a coleta cheia (n=19200) com o mesmo
# mecanismo (scripts/collect_traces_manual.sh).
#
# Captura, por runtime x escalonador x algoritmo (8 runs, 1 rep):
#   StarPU (dmda,dmdas): FxT (prof_file_*)            + STARPU_BUS_STATS/WORKER_STATS no log
#   PaRSEC (lfq,gd)    : .prof (PINS) + DAG .dot       + device_show_statistics  no log
# ou seja: tudo que o gpu_doe_sweep.sh ja pegava + os CONTADORES NATIVOS de
# transferencia H2D/D2H (TRACE_STATS=1) -- a evidencia direta do thrashing.
#
#   uso (dentro de uma alocacao interativa no poti):
#     bash scripts/pcad/collect_validate.sh
#
#   knobs (env):  N (default 7680)  B (960)  THREADS (16)  GPUS (1)
#                 NIX_CMD ("nixw nix" na PCAD; use "nix" localmente)
#
# Saida: data/validate_<stamp>/{results.csv, design.csv, runs/<...>/}
# Depois: sincronize de volta (scripts/pcad/get_from_pcad.sh) e rode o
# pos-processamento + plots listados no fim deste script.

set -euo pipefail
cd "$(dirname "$0")/../.."   # raiz do repo

# Caso pequeno: rapido, mas grande o bastante (8x8 tiles em b=960) para a GPU
# offloadar gemms e os padroes de transferencia/dependencia aparecerem.
N="${N:-7680}"
B="${B:-960}"
THREADS="${THREADS:-16}"
GPUS="${GPUS:-1}"
NIX_CMD="${NIX_CMD:-nixw nix}"

STAMP="$(date +%Y%m%d_%H%M%S)"
RESULTS_DIR="${RESULTS_DIR:-data/validate_${STAMP}}"
# NB: nome != design.csv -- collect_traces_manual.sh copia o DoE para
# $RESULTS_DIR/design.csv, entao o fonte tem que ser outro arquivo.
DESIGN_FILE="${RESULTS_DIR}/validate_doe.csv"

mkdir -p "$RESULTS_DIR"

# DoE minimo: 2 algoritmos x 2 runtimes x seus 2 escalonadores, FP32, 1 rep.
cat > "$DESIGN_FILE" <<EOF
precision,algorithm,n,b,runtime,scheduler,rep
FP32,potrf,${N},${B},starpu,dmda,1
FP32,potrf,${N},${B},starpu,dmdas,1
FP32,geqrf,${N},${B},starpu,dmda,1
FP32,geqrf,${N},${B},starpu,dmdas,1
FP32,potrf,${N},${B},parsec,lfq,1
FP32,potrf,${N},${B},parsec,gd,1
FP32,geqrf,${N},${B},parsec,lfq,1
FP32,geqrf,${N},${B},parsec,gd,1
EOF

echo "==> Validacao: N=$N B=$B t=$THREADS g=$GPUS"
echo "==> Saida    : $RESULTS_DIR"
echo "==> DoE      : $DESIGN_FILE (8 runs)"
echo

# Reusa o orquestrador validado: entra nos dois dev shells e roda o sweep com
# TRACE/TRACE_FULL/TRACE_DAG/TRACE_STATS ligados.
RESULTS_DIR="$RESULTS_DIR" DESIGN_FILE="$DESIGN_FILE" \
  THREADS="$THREADS" GPUS="$GPUS" NIX_CMD="$NIX_CMD" TRACE_STATS=1 \
  bash scripts/collect_traces_manual.sh both

echo
echo "=================================================================="
echo " Coleta de validacao pronta em: $RESULTS_DIR"
echo
echo " Checagem rapida (deve existir, por run):"
echo "   StarPU: prof_file_*               + 'Data transfer stats:' no .log (NUMA<->CUDA = H2D/D2H)"
echo "   PaRSEC: cham_<k>-0.prof-* + .dot  + 'All Devs' no .log (pares Required/Transfered = H2D/D2H)"
echo "   grep -l 'Data transfer stats'  $RESULTS_DIR/runs/*starpu*/*.log"
echo "   grep -l 'All Devs'             $RESULTS_DIR/runs/*parsec*/*.log"
echo
echo " Proximos passos (no shell de analise, apos sincronizar de volta):"
echo "   # 1) StarPU FxT -> parquets do StarVZ"
echo "   nix develop --command bash scripts/analysis/trace_phase1.sh $RESULTS_DIR/runs/*starpu*"
echo "   # 2) PaRSEC .prof+.dot -> tasks/dag parquets"
echo "   export DBP2XML=\$(nix build --impure --no-link --print-out-paths .#parsec)/bin/dbp2xml"
echo "   nix develop --command bash -c '"
echo "     Rscript scripts/analysis/parsec_tasks_to_parquet.r $RESULTS_DIR/runs/*parsec*"
echo "     Rscript scripts/analysis/parsec_dag_to_parquet.r   $RESULTS_DIR/runs/*parsec*'"
echo "   # 3) Plots (passando o base_dir desta coleta)"
echo "   for s in plot_task_times plot_bounds plot_wait_time plot_scheduler_health; do"
echo "     nix develop --command Rscript scripts/analysis/\$s.r $RESULTS_DIR/runs; done"
echo "=================================================================="
