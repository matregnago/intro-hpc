# Plano de experimentos — impacto do escalonador na GPU (StarPU vs PaRSEC)

Estudo do **impacto dos escalonadores** em Cholesky (`potrf`) e QR (`geqrf`) na
GPU, comparando os runtimes StarPU e PaRSEC sob o Chameleon. É a adaptação,
para um **nó único com GPU**, da metodologia de Schnorr/Marcondes
(WAMCA 2025, *"Impact of Data Distribution and Schedulers for the LU
Factorization on Multi-Core Clusters"* — ver `contexto/wamca-2025/`).

## Como mapeia no trabalho do Schnorr

| Schnorr (multi-nó CPU, StarPU-MPI) | Aqui (nó único, GPU) |
|---|---|
| transporte MPI {nmad, openmpi} | **runtime {StarPU, PaRSEC}** |
| scheduler {lws,random,dmda,dmdas} | **StarPU {dmda, dmdas}** vs **PaRSEC {lfq, gd}** |
| distribuição de dados PxQ (só com MPI) | **precisão {FP32, FP64}** (nossa fonte de heterogeneidade) |
| LU (`dgetrf_nopiv`) | **Cholesky (`potrf`) + QR (`geqrf`)** |
| volume de comunicação MPI (métrica causal) | **tráfego H2D/D2H na GPU** (thrashing do pushout do DTD do PaRSEC) |

Decisões de escopo (travadas): validar **local primeiro** (RTX 4060 Ti, WSL2) e
levar para o **PCAD poti** (RTX 4070) depois; 2º fator = **FP32 vs FP64**;
escalonadores **dmda/dmdas (StarPU) vs lfq/gd (PaRSEC)**.

## Etapas

### E1 — Seleção de tile (`b`)
Análogo ao *block-size* do Schnorr. Fixa N e o escalonador campeão de cada
runtime (dmda / lfq), varia só o tile, por (runtime × algoritmo × precisão).
Ancora o `b` usado em E2 e já justifica a escolha (crítica recorrente de
revisor).

```bash
nix develop --impure --command Rscript scripts/doe/doe_gpu_tile.r   # -> scripts/doe/doe_gpu_tile.csv
DESIGN_FILE=scripts/doe/doe_gpu_tile.csv bash scripts/run_gpu_sched.sh
```

### E2 — Varredura de escalonadores (experimento principal)
Fixa N e `b` (de E1), cruza **scheduler × runtime × algoritmo × precisão**, 10
replicações randomizadas (DoE.base). É a figura-mãe.

```bash
nix develop --impure --command Rscript scripts/doe/doe_gpu_sched.r  # -> scripts/doe/doe_gpu_sched.csv
bash scripts/run_gpu_sched.sh
# saída: data/gpu_sched_<timestamp>/{results.csv,design.csv,runs/,plots/}
```

Knobs por env: `N`, `B`, `REPS` (geração do DoE); `THREADS`, `GPUS`, `DESIGN_FILE` (execução).

### E3 — Escalabilidade em N
Melhor escalonador de cada runtime, varrendo N. Usa a infra que já existe:
`scripts/run_local_gpu_compare.sh` + `scripts/analysis/plot_local_gpu_compare.r`
(GFLOPS×N + razão PaRSEC/StarPU).

### E4 — Rastros (StarVZ) para explicar o *porquê*
Configs representativas (melhor vs pior escalonador; StarPU-dmda vs PaRSEC-lfq).
Painéis Espaço-Tempo **com as lanes de GPU** (já habilitadas no PaRSEC). Métricas
causais: contagem/tempo de H2D-D2H, idle por worker, gaps, duração média por
tipo de tarefa, % de tarefas na GPU vs CPU. É onde o **thrashing de pushout do
DTD do PaRSEC** (ver `docs/cuda_parsec.md`) vira resultado. Reaproveita
`scripts/analysis/plot_traces.r` (StarVZ lê StarPU e PaRSEC).

## Como mostrar os dados

1. **Violino de tile** (E1) — makespan × `b`, facet por algoritmo.
2. **Violino principal** (E2) — `makespan_by_scheduler.{png,pdf}`: makespan ×
   escalonador, cor = runtime, `facet_grid(algoritmo ~ precisão)`. Cópia do
   `distrib-scheduler.pdf` do Schnorr. Versão em GFLOPS junto.
3. **Linha GFLOPS × N + barras de razão** (E3) — já existe.
4. **Painéis Espaço-Tempo** (E4) — StarVZ `panel_st` com lanes de GPU.
5. **Tabelas** — `summary_gpu_sched.csv` (média ± dp, CV) e
   `best_scheduler.csv` (campeão por célula); métricas de rastro em E4.

## Como analisar as diferenças

- **DoE formal** (`DoE.base::fac.design`, randomizado, replicado) → rigor.
- **Estatística**: média ± dp e CV por célula; ANOVA do makespan com fatores
  `scheduler*runtime*algoritmo*precisão`, ou comparações pareadas; violinos com
  dispersão visível.
- **Explicação causal**: correlacionar o ranking de makespan com as métricas de
  rastro (transfers H2D/D2H, idle) — não só *quem* ganha, mas *por quê*.

## Arquivos desta infra

- `scripts/doe/doe_gpu_sched.r` / `doe_gpu_tile.r` / `doe_gpu_traces.r` — geram os DoE (E2 / E1 / E4).
- `scripts/gpu_doe_sweep.sh` — sweep DoE-driven, por runtime, dentro do dev shell; honra `TRACE=1` (FxT StarPU / profile PaRSEC) para a E4.
- `scripts/run_gpu_sched.sh` — orquestrador local (entra nos dois shells + plota).
- `scripts/analysis/plot_gpu_sched.r` — violinos + tabelas-resumo (E2).
- `slurm/gpu_phases.slurm` — job no PCAD (`poti`) que roda E1→E4 no molde do `traces.slurm` (rsync→`$SCRATCH`, `nixw nix develop`, trap de copyback), com `TRACE=1` na E4.

## Rodar no PCAD

```bash
# 1) gerar os DoE localmente (shell default, tem R + DoE.base)
Rscript scripts/doe/doe_gpu_tile.r     # E1 -> scripts/doe/doe_gpu_tile.csv
Rscript scripts/doe/doe_gpu_sched.r    # E2 -> scripts/doe/doe_gpu_sched.csv
Rscript scripts/doe/doe_gpu_traces.r   # E4 -> scripts/doe/doe_gpu_traces.csv
# 2) sincronizar e submeter
bash scripts/pcad/sync_to_pcad.sh
ssh pcad 'cd ~/intro_hpc && sbatch slurm/gpu_phases.slurm'
# 3) trazer de volta + plotar
bash scripts/pcad/get_from_pcad.sh
# resultados em data/chameleon-gpu-phases_<jobid>/{e1_tile,e2_sched,e3_scal,e4_traces}/
```

Os rastros da E4 ficam em `data/.../e4_traces/runs/<id>/` (FxT do StarPU + `cham_*` do
PaRSEC), prontos para o `scripts/analysis/plot_traces.r` (StarVZ).

> No PCAD, todo `nix` precisa do shim `nixw` (ver `CLAUDE.md`); as fases do job
> são independentes (uma falha não aborta as demais).


k - iteration plot
tempo de cada tarefa
validar se o resultado do chameleon parsec funciona
