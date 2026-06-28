# Painéis StarVZ para comparar StarPU × PaRSEC: viabilidade e scripts

Análise das ideias do `graficos.md` (recomendadas pelo Schnorr) — o que é mensurável em
**cada runtime** dado o que o PaRSEC expõe — e os scripts de plotagem implementados. Os
requisitos de cada painel foram conferidos no fonte do StarVZ v0.8.4 (`phase2_kchart.R`,
`phase2_states_chart.R` para ABE/CPB, `phase1-workflow.sh` para o caminho PaRSEC). Contexto de
como o PaRSEC vira os tibbles do StarVZ: [`starvz_e_captura_parsec.md`](starvz_e_captura_parsec.md).

## Viabilidade das 4 ideias

| Ideia | StarPU | PaRSEC | Observação |
|---|:---:|:---:|---|
| **1. Tempo médio por tipo de tarefa** | ✅ | ✅ hoje | só precisa de `application.parquet` (kernel + duração) |
| **2. Bound ABE** | ✅ | ✅¹ | `geom_abe` usa só durações (não usa DAG/GFlop/pmtool) |
| **3. K-iteration** | ✅ | ❌ | exige `Application$Iteration`; no PaRSEC é `NA` (não há índice de iteração no trace) |
| **4. Wait time da tarefa** | ✅ | ✅² | `StartTime − max(EndTime dos predecessores)`, do DAG |

¹ **ABE no PaRSEC tem uma pegadinha:** o painel embutido do StarVZ (`abe_cpu_cuda` dentro do
`panel_st`) **quebra** no PaRSEC porque depende do `ResourceType` ∈ {CPU, CUDA}, e o
`dbp2paje` corrompe isso para `"MVT-"`. Por isso o ABE *no gantt* (`plot_traces.r`) fica só no
StarPU, mas o ABE **calculado direto das durações** (`plot_bounds.r`) funciona nos dois — é a
mesma área `Σ duração / nº workers`.

² **Wait time no PaRSEC** precisa do `.dot` do grapher + `tasks.parquet` reconstruído
(`parsec_dag_to_parquet.r` / `parsec_tasks_to_parquet.r`). Caveat honesto: é *prontidão por
dados* (elegibilidade pelo DAG), **não** a latência exata da fila do escalonador.

**K-iteration (ideia 3) ficou de fora.** Para suportá-lo no PaRSEC seria preciso *derivar* o
índice de iteração da estrutura do DAG (camada topológica) — trabalho extra e aproximado.

## Scripts implementados

Todos em `scripts/analysis/`, reusam `trace_common.r` (parse do nome do run, normalização de
kernel `dgemm`↔`gemm`, descoberta de runs, **fonte unificada** `read_exec()` e normalização de
tempo p/ ms). Default = **job de GPU** `data/manual_traces_20260621_172820/runs` (FP32
spotrf/sgeqrf n19200 b960; 2 runtimes × 2 escalonadores × 2 algoritmos). Saídas em `plots/`
(PNG 300dpi + PDF + CSV). Todos facetam por algoritmo e aceitam `base_dir` como 1º argumento.

| Script | Ideia | Fonte | GPU? |
|---|---|---|---|
| `plot_task_times.r` | 1 + violino | `read_exec` (application StarPU / tasks PaRSEC) | ✅ os dois |
| `plot_bounds.r` | 2 + add-4 | `read_exec` (makespan/ABE) + `dag.parquet` (CPB) | ✅ makespan+CPB; ABE só CPU |
| `plot_wait_time.r` | 4 | `tasks.parquet` (normalizado p/ ms) | ✅ os dois |
| `plot_scheduler_health.r` | add-2/3 | `read_exec` (concorrência) + `variable`/`tasks` (Ready) | ✅ os dois |
| `plot_transfers.r` | thrashing H2D/D2H | contadores nativos no `.log` (`TRACE_STATS=1`) | ✅ os dois |
| `plot_traces.r` (estendido) | 2 (visual) | `starvz_read` | gantt StarPU (PaRSEC GPU não carrega) |

> `plot_transfers.r` precisa da coleta nova (com `TRACE_STATS=1`); os logs antigos
> (`manual_traces`) não têm os contadores.

### Pré-processamento do job de GPU (uma vez)

O job de GPU vem com trace **bruto**; converter antes de plotar:

```bash
# StarPU: FxT -> parquets do StarVZ (application/dag/tasks/variable)
nix develop --command bash scripts/analysis/trace_phase1.sh \
  data/manual_traces_20260621_172820/runs/000[1-4]_starpu_*
# PaRSEC: .prof + .dot (grapher) -> tasks/dag.parquet (precisa de dbp2xml no PATH;
#   ele vem do `nix build --impure .#parsec`)
export DBP2XML=$(nix build --impure --no-link --print-out-paths .#parsec)/bin/dbp2xml
nix develop --command bash -c '
  Rscript scripts/analysis/parsec_tasks_to_parquet.r data/manual_traces_*/runs/000[5-8]_parsec_*
  Rscript scripts/analysis/parsec_dag_to_parquet.r   data/manual_traces_*/runs/000[5-8]_parsec_*'
```

### Plotar (dentro do `nix develop`)

```bash
Rscript scripts/analysis/plot_task_times.r        # tempo médio + violino por kernel
Rscript scripts/analysis/plot_bounds.r            # makespan vs CPB (e ABE)
Rscript scripts/analysis/plot_wait_time.r         # wait time por tarefa
Rscript scripts/analysis/plot_scheduler_health.r  # concorrência + Ready/starvation
# Para o job CPU (dpotrf n19200 b480), passe o base_dir:
Rscript scripts/analysis/plot_task_times.r data/parcial/chameleon-runtimes_782088/runs
```

## Coleta da rodada nova (com contadores nativos de transferência)

A coleta atual (`gpu_doe_sweep.sh` via `collect_traces_manual.sh`) já pega FxT (StarPU),
`.prof`+`.dot` (PaRSEC) e — com o knob **`TRACE_STATS=1`** (novo) — grava os **contadores
nativos H2D/D2H no `.log`** de cada run: StarPU `STARPU_BUS_STATS` (bloco "Data transfer
stats:" NUMA↔CUDA) e PaRSEC `device_show_statistics` (linha "All Devs" com pares
Required/Transfered). É a evidência direta do thrashing, e é o que o `plot_transfers.r` lê.

**Validação interativa no PCAD** (poti) antes da coleta cheia:

```bash
bash scripts/pcad/collect_validate.sh      # N pequeno (7680), 8 runs, TRACE_STATS=1
# checagem: grep -l 'Data transfer stats' .../runs/*starpu*/*.log
#           grep -l 'All Devs'            .../runs/*parsec*/*.log
```

Validado localmente ponta-a-ponta (RTX 4060 Ti, N=2880): coleta → logs com contadores →
phase1 + adaptadores + os 5 plots, tudo OK. O D2H/kernel do PaRSEC ≈ um tile inteiro (3,5 MB,
writeback completo) vs StarPU bem menor (residência) — já visível mesmo no caso pequeno. A
coleta cheia usa o mesmo mecanismo com N=19200 (`collect_traces_manual.sh`).

## Pegadinhas importantes (GPU)

1. **Unidades.** Os parquets do StarVZ (StarPU, e PaRSEC via dbp2paje) estão em **ms**; o
   adaptador do PaRSEC (`parsec_*_to_parquet.r`) emite **µs**. `trace_common.r::time_scale_to_ms`
   normaliza tudo p/ ms (heurística: tem `application.parquet` ⇒ ms; só `tasks.parquet` ⇒ µs).
   Sem isso a comparação fica off por 1000×. Validação: kernels de QR ficam comparáveis entre
   runtimes (ambos ~ms); o gap de 1000× sumiu.
2. **ABE não vale na GPU.** O ABE (limite por área = trabalho/Pworkers) é um conceito
   *homogêneo*. Na GPU o PaRSEC offloada pras streams mas o `tasks.parquet` atribui o trabalho a
   ~16 workers CPU → o ABE estoura o makespan (eff < 1, impossível p/ limite inferior). Use o
   **CPB (caminho crítico)**, que é válido nos dois (todos eff_cpb > 1). ABE fica como
   diagnóstico só no caso CPU homogêneo.
3. **Fonte por runtime na GPU.** StarPU usa `application.parquet` (split CPU/CUDA limpo); o
   PaRSEC GPU **não** tem `application` (o `pj_dump` do dbp2paje falha no paje de GPU) e usa
   `tasks.parquet`. As lanes de stream da GPU **não** entram no `states.parquet` (o parser só
   casa `Thread N`, não `Stream N`) — mas o `tasks.parquet` é GPU-completo (timing via
   `tpid:tid` sobre todos os workers), então é a fonte certa.
4. **Bug corrigido:** `parsec_dag_to_parquet.r` não limpava `LD_LIBRARY_PATH` ao chamar
   `dbp2xml` (o `parsec_tasks_to_parquet.r` já fazia) → `dbp2xml` falhava só no conversor de DAG.

## Resultados de GPU (sanity)

- **spotrf FP32**: StarPU ~250 ms vs PaRSEC ~1040 ms (~4×) — gap dominado pelo thrashing D2H do
  PaRSEC (gemm ~15 ms/tarefa no PaRSEC vs ~0.2 ms no StarPU; ver `docs/cuda_parsec.md`). **sgeqrf**
  fica mais perto (StarPU ~2.4 s vs PaRSEC ~3.6 s) pq o painel é na CPU.
- **CPB** (eff makespan/CPB): spotrf StarPU ~2.2–2.5, PaRSEC ~1.7–1.9; sgeqrf ~2.0–2.5 (StarPU),
  ~1.2–1.4 (PaRSEC).
- **Starvation** (Ready==0): spotrf StarPU ~62–65 % (GPU termina os tiles rápido e a fila CPU
  esvazia) vs PaRSEC ~14–22 %.

## Limitação restante

Os runs de GPU têm 1 réplica e o split CPU/GPU do PaRSEC se perde no adaptador (tudo cai em
workers CPU). Para split de dispositivo confiável no PaRSEC e mais réplicas, re-coletar no
PCAD (poti/RTX 4070). O job `data/parcial` (dpotrf n19200 b480) é **CPU-only** e serve de
contraponto homogêneo (onde o ABE é válido, eff ~1.04).
