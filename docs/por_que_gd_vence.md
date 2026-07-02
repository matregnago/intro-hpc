# Por que o `parsec:gd` aparece como melhor escalonador no poti (e por que isso não contradiz nada)

Documento de evidências para a figura *GFLOPS vs N* do poti (FP64, b=500), onde
`parsec:gd` lidera sobre `starpu:dmda`/`dmdas` apesar de (i) o PaRSEC não ter
escalonador GPU-aware (decisão CPU-vs-GPU estática no device layer — ver
[`analise_escalonadores_gpu.md`](analise_escalonadores_gpu.md)) e (ii) o StarPU
ter perf-model medido por worker. Todos os números vêm de
`data/gpu_traces_poti_20260701_131857` (dpotrf/dgetrf_nopiv, N=40000, b=500),
`data/gpu_{n_size,tile}_{poti,tupi}_*` e das tabelas em
`plots/final/why_gd_poti/`.

## TL;DR

**O `gd` não vence porque lida bem com GPU; vence porque no poti+FP64 a jogada
certa é *usar menos* a GPU — e o PaRSEC faz isso por acidente de design,
enquanto o StarPU paga um custo estrutural de CPU para saturá-la.** Nos regimes
em que a GPU vale a pena (FP32 no poti; FP64 na RTX 4090 do tupi), o quadro
inverte e o `dmda`/`dmdas` vence com folga — até 2.3×.

## 1. O regime de hardware: em FP64 a CPU do poti é o recurso majoritário

A RTX 4070 é Ada de consumo: FP64 castrado em 1/64 do FP32 → **pico FP64 =
458 GFlops** (o PaRSEC imprime isso no log; ver
[`parsec_cuda_results.md`](parsec_cuda_results.md)). Capacidades **efetivas
medidas** nos traces (`gpu_split_summary.csv`, dpotrf):

| classe | starpu:dmda | parsec:gd |
|---|---:|---:|
| GPU (GFlops com o device ativo) | 415 | 392 |
| CPU (19 workers × taxa por worker ocupado) | **605** | **763** |

Em double, as 19 CPUs valem ~1.5–1.9 GPUs. Em FP32 a mesma GPU tem 29.3 TFlops
(~30× a CPU) — regime oposto.

## 2. A GPU é o gargalo em todos os runs (bounds)

`plots/final/gpu_traces_poti_20260701_131857/bounds_summary.csv`: `eff_gpu`
(makespan / tempo-de-device-ativo) fica em **1.01–1.02 para os 4
escalonadores**, nos dois algoritmos. A GPU trabalha de ponta a ponta e define
o makespan; a corrida se decide por *quanto* cada runtime carrega a GPU e
*quanto de CPU sobra para o resto*.

## 3. O split de trabalho: dmda maximiza a GPU, gd maximiza a CPU

Figuras `gpu_split_contrib` e `gpu_split_share` (tabela
`gpu_split_summary.csv`). GFlops entregues por classe (gflops da classe /
makespan), dpotrf N=40000:

| config | CPU | GPU | total | % flops na GPU | r* do próprio run |
|---|---:|---:|---:|---:|---:|
| parsec:gd | **663** | 392 | **1054** | 37.1% | 34.0% |
| starpu:dmdas | 591 | 416 | 1007 | 41.3% | 41.1% |
| starpu:dmda | 593 | 412 | 1005 | 41.0% | 40.7% |
| parsec:lfq | 553 | 391 | 944 | 41.4% | 34.2% |

(dgetrf_nopiv segue o mesmo padrão: gd 656+406=1062 vs dmda 595+416=1011.)

Duas leituras importantes:

- **O dmda não está "errando" o balanceamento**: ele opera exatamente no split
  ótimo *do sistema que ele enxerga* (41.0% vs r*=40.7%). O perf-model faz o
  trabalho dele.
- **O problema é que o sistema do StarPU é pior**: capacidade total dmda =
  605+415 = **1020 GFlops** — *abaixo* dos 1054 que o gd entrega. Mesmo com
  escalonamento perfeito o dmda não alcança o gd. A perda é estrutural, não de
  política de escalonamento.

## 4. De onde vem a perda estrutural do StarPU (duas fontes)

### 4a. CPU 26% mais lenta por tarefa (task_times_summary.csv, classe CPU)

| gemm b=500 em CPU | média | mediana |
|---|---:|---:|
| starpu:dmda | 7.81 ms | 7.36 ms |
| parsec:gd | **6.19 ms** | **4.09 ms** |

Mecanismo: o StarPU dedica **4 threads de worker CUDA**
(`STARPU_NWORKER_PER_CUDA=4`) + thread principal além dos 19 workers de CPU →
24 threads em 20 cores (8P+12E). O PaRSEC gerencia a GPU
*oportunisticamente* (o worker que pega o mutex do device progride os streams,
sem thread dedicada — `dev_cuda.c:1991`), então seus 19 workers rendem mais
por segundo ocupado (39.7 vs 31.9 GFlops/worker·s). Ver figura
`task_times_mean`/`task_times_violin` (agora particionadas por classe CPU/GPU).

### 4b. Thrashing de VRAM no caminho H2D (transfers_summary.csv)

Matriz N=40000 FP64 = **12.8 GB**; VRAM útil = **10.4 GiB** — a matriz não
cabe. Contadores nativos (dpotrf):

| config | H2D | D2H |
|---|---:|---:|
| starpu:dmda | **36.4 GB** (~2.8× a matriz) | 31.1 GB |
| starpu:dmdas | 32.5 GB | 26.7 GB |
| parsec:gd | **6.0 GB** (3.4% do requerido — cache LRU) | **61.4 GB** (100% — pushout) |

O StarPU re-encena a matriz ~3× por evicção; esse H2D fica no caminho crítico
dos kernels (o kernel espera o input). O PaRSEC quase não re-encena (LRU de
device segura 96.6% dos inputs), mas paga writeback eager de **todo** tile
escrito (o `pushout` incondicional da nossa camada DTD — 1.86 MB/kernel,
`cuda_parsec.md`); D2H, porém, sai do caminho crítico (stream dedicado, o
kernel seguinte não espera). Em dgetrf_nopiv: StarPU 76.3/66.3 GB vs PaRSEC
12.0/123.8 GB.

## 5. Por que `gd` > `lfq` dentro do PaRSEC

- `lfq` sobrecarrega mais a GPU (41.4% dos flops vs 37.1%; 36539 vs 32964
  kernels no log) — e a GPU é o gargalo (§2);
- `lfq` tem variância bimodal (CV ~7.7% vs ~2.6% do gd): work-stealing
  reordena a chegada das tarefas ao hook CUDA e muda o tile residente
  "vencedor" a cada run — ver
  [`variancia_parsec_gpu.md`](variancia_parsec_gpu.md) §2a/§3.

## 6. As contraprovas: onde a GPU vale a pena, o StarPU esmaga

Figura `precision_grid` (tabela `precision_grid_summary.csv`) — melhor média
do sweep (n,b) por config, três regimes ordenados por pico de GPU:

| regime | melhor PaRSEC | melhor StarPU | quem lidera |
|---|---:|---:|---|
| poti RTX 4070 · FP64 (pico GPU 458 GF) | gd **1075** (getrf) / **1055** (potrf) | dmdas 998 / 991 | **parsec:gd** |
| tupi RTX 4090 · FP64 (pico ~1.3 TF) | gd 1611 / 1466 | dmdas **1837** / dmda **1758** | **starpu** |
| poti RTX 4070 · FP32 (pico 29.3 TF) | lfq 6432 / 7949 | dmdas **14929** / **14689** | **starpu, ~2×** |

O `gd` só lidera no único regime em que o pico FP64 da GPU fica **abaixo** do
agregado das CPUs. Assim que a GPU passa a valer a pena (4090 em FP64, ou
FP32), o perf-model + prefetch do `dmda` fazem diferença real e o device layer
estático do PaRSEC (peso teórico invertido + contador de carga, sem tempo
medido — `analise_escalonadores_gpu.md` §4) fica para trás. No próprio poti
FP64, o StarPU ainda vence onde a pressão de VRAM some (N≤10000; b=250).

## 7. A frase para a apresentação

> No poti, em dupla precisão, a RTX 4070 entrega no máximo 458 GFlops — menos
> do que as 19 CPUs. Nesse regime o escalonador ótimo usa a GPU como acessório,
> e é isso que o PaRSEC faz *por construção* (heurística estática que
> under-offloada) e o `gd` faz melhor ainda (CPU-throughput máximo, sem threads
> dedicadas de GPU). O StarPU `dmda` balanceia perfeitamente o sistema que
> enxerga, mas o sistema dele é pior: 4 threads roubadas para os workers CUDA
> (-26% por tarefa de CPU) e 3 matrizes re-encenadas por evicção. Quando a GPU
> passa a valer a pena — FP32, ou a 4090 do tupi — o quadro inverte e o
> `dmda` vence por até 2.3×. `gd` não é "melhor escalonador de GPU"; é o
> escalonador que menos atrapalha a CPU num nó cuja GPU é fraca em FP64.

## Figuras/tabelas (em `plots/final/why_gd_poti/`)

| arquivo | mostra |
|---|---|
| `gpu_split_contrib.{png,pdf}` | contribuição CPU vs GPU nos GFlops totais (o plot-chave) |
| `gpu_split_share.{png,pdf}` | % dos flops na GPU vs split ótimo r* por run |
| `gpu_split_rates.{png,pdf}` | capacidade efetiva por classe (CPU StarPU 26% abaixo) |
| `precision_grid.{png,pdf}` | os três regimes: gd só lidera em poti+FP64 |
| `transfers_d2h.{png,pdf}` | H2D/D2H vs tamanho da matriz e VRAM (thrashing vs pushout) |
| `task_times_{mean,violin}.{png,pdf}` | duração por kernel × classe (gemm CPU 7.8 vs 6.2 ms) |
| `*_summary.csv` | tabelas por trás de cada figura |

Scripts: `scripts/analysis/plot_gpu_split.r`, `plot_precision_grid.r` (novos),
`plot_transfers.r` (parse PaRSEC corrigido p/ unidades B/MB/GB),
`plot_task_times.r` (split por classe de recurso). Reproduzir com
`PLOTS_DIR=plots/final/why_gd_poti` apontando para
`data/gpu_traces_poti_20260701_131857/runs`.
