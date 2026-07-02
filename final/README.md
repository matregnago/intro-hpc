# final/ — pipeline da versão final (estudo de escalonador na GPU)

Bundle auto-contido dos scripts usados para gerar os resultados finais do estudo
StarPU vs PaRSEC na GPU (Chameleon, nó único com GPU). São cópias congeladas das
versões validadas — os originais continuam em `scripts/` para os fluxos antigos.

## Conteúdo

| Arquivo | Papel |
|---|---|
| `doe_gpu_tile.r` | Gera os DoEs de **E1 (block size)**, um por nó → `doe_gpu_tile_{poti,tupi}.csv` |
| `doe_gpu_tile_poti.csv` / `_tupi.csv` | DoE de E1 (~360 runs cada; ver parâmetros) |
| `doe_gpu_sched.r` | Gera o DoE de **E2 (impacto do escalonador)** → `doe_gpu_sched.csv` |
| `doe_gpu_trace_ideal.r` | Gera o DoE de **E4 (rastro StarVZ no tile ideal)** → `doe_gpu_trace_ideal_{poti,tupi}.csv` |
| `gpu_doe_sweep.sh` | Motor do sweep (calibra dmda/dmdas, roda as runs, grava `time,gflops`) |
| `run_gpu_sched.sh` | Orquestrador **local** (roda os 2 runtimes + plota) |
| `gpu_tile_{poti,tupi}.slurm` | Job SLURM de **E1** (poti=RTX 4070 THREADS=19; tupi=RTX 4090 THREADS=23), `TRACE=0` |
| `gpu_trace_ideal_{poti,tupi}.slurm` | Job SLURM de **E4** (mesmo motor, `TRACE=1 TRACE_FULL=1 TRACE_DAG=1 TRACE_STATS=1`) |

> **GPU workers (E1 e E4).** O `gpu_doe_sweep.sh` roda o StarPU com
> `STARPU_NWORKER_PER_CUDA=4` (default, sobrescrevível por env) para casar com os
> 4 exec streams do PaRSEC. Vale para a **calibração** e para as runs cronometradas,
> nos dois jobs — assim o tile ideal de E1 é medido na mesma config do rastro de E4
> e os painéis Espaço-Tempo do StarVZ têm o mesmo nº de lanes de GPU nos dois runtimes.

## Parâmetros fixados — E1 (block size)

O objetivo de E1 é traçar a **montanha** de GFLOPS vs `b` e ler o **pico** = block
size ideal (que ancora `b` em E2/E4). Tamanhos grandes o bastante para **saturar
a GPU** (o design antigo N=20000/N=40000 deixava o FP32 curto demais e o pico na
borda da grade — ver histórico no cabeçalho de `doe_gpu_tile.r`):

- **N por precisão** (a matriz inteira fica no host; poti tem 96 GB DDR5, RAM
  nunca é gargalo):
  - **FP32 → N = 80000**, `b ∈ {1000, 2000, 3200, 4000, 5000, 8000}` (N/b 80→10;
    denso em 3200-5000, onde o pico deve cair).
  - **FP64 → N = 40000**, `b ∈ {250, 320, 400, 500, 625, 800, 1000, 2000, 4000}`
    (N/b 160→10; adensado **abaixo de 1000** — a varredura antiga `{1000,2000,4000}`
    achava o pico sempre na borda `b=1000`, logo o topo da montanha FP64 caía fora
    da grade).
- `b` **sempre divisor de N** — PaRSEC+CUDA crasha (SIGSEGV/heap-corrupt) em tile
  não-divisor (job 797989).
- **Fatorações**: `potrf` (Cholesky) + `getrf_nopiv` (LU). QR está fora (resultado
  numérico inválido no PaRSEC+GPU).
- **Escalonadores**: StarPU `{dmda, dmdas}`, PaRSEC `{lfq, gd}`.
- **THREADS** por nó (reserva 1 core p/ a GPU): poti=19, tupi=23.
- **reps = 3**, ordem das linhas randomizada. E1 **não coleta trace** — só
  `time,gflops` (o que ranqueia o tile).
- Binding de CPU **ligado** (Linux nativo); desligado automaticamente só sob WSL.

`N` e tiles são iguais nos dois nós (comparação 4070 vs 4090 no mesmo problema);
o que muda por nó é só o `THREADS`. ~360 runs/nó (FP32 144 + FP64 216);
`--time=06:00:00`.

## Como rodar no PCAD

```bash
# 1. (do PC local) sincronize o repo
bash scripts/pcad/sync_to_pcad.sh

# 2. (no front-end gppd-hpc) regenere os CSVs se preciso
nix develop --impure --command Rscript final/doe_gpu_tile.r

# 3. E1 — submeta um job de block size por partição
sbatch final/gpu_tile_poti.slurm
sbatch final/gpu_tile_tupi.slurm

# 4. (do PC local) puxe os resultados e leia o pico
bash scripts/pcad/get_from_pcad.sh
Rscript scripts/analysis/plot_block_size.r data/gpu_tile_poti_<jobid>
#   -> plots/gflops_vs_b_r.{png,pdf} + plots/block_size_peak.csv (b ideal)
```

Resultados em `data/gpu_tile_{poti,tupi}_<jobid>/results.csv` (+ modelo StarPU
calibrado em `.starpu/`). Depois de E1: escolha o `b` campeão por kernel, gere E2
com `final/doe_gpu_sched.r` e, para os rastros, atualize `doe_gpu_trace_ideal.r`
com o `b` ideal novo antes de submeter `final/gpu_trace_ideal_{poti,tupi}.slurm`.

## Smoke test local

```bash
# N=80000 é pesado localmente; reduza por env (ex.: gere um CSV com N menor) ou
# rode poucas linhas. O caminho em si:
DESIGN_FILE=final/doe_gpu_tile_poti.csv THREADS=8 \
  RESULTS_DIR=data/gpu_tile_local bash final/run_gpu_sched.sh
```
