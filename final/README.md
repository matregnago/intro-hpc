# final/ — pipeline da versão final (estudo de escalonador na GPU)

Bundle auto-contido dos scripts usados para gerar os resultados finais do estudo
StarPU vs PaRSEC na GPU (Chameleon, nó único com GPU). São cópias congeladas das
versões validadas — os originais continuam em `scripts/` para os fluxos antigos.

## Conteúdo

| Arquivo | Papel |
|---|---|
| `doe_gpu_tile.r` | Gera os DoEs de **E1 (block size)**, um por nó → `doe_gpu_tile_{poti,tupi}.csv` |
| `doe_gpu_tile_poti.csv` / `_tupi.csv` | DoE de E1 (240 runs cada: 2 prec × 2 fatorações × 5 tiles × 4 sched × 3 reps) |
| `doe_gpu_sched.r` | Gera o DoE de **E2 (impacto do escalonador)** → `doe_gpu_sched.csv` |
| `gpu_doe_sweep.sh` | Motor do sweep (calibra dmda/dmdas, roda as runs, grava `time,gflops`) |
| `run_gpu_sched.sh` | Orquestrador **local** (roda os 2 runtimes + plota) |
| `gpu_tile_poti.slurm` | Job SLURM de E1 na partição **poti** (RTX 4070, THREADS=19) |
| `gpu_tile_tupi.slurm` | Job SLURM de E1 na partição **tupi** (RTX 4090, THREADS=23) |

## Parâmetros fixados

- **N = 20000** nos dois nós (comparação 4070 vs 4090 no mesmo problema).
- **Tiles**: `b ∈ {400, 500, 800, 1000, 2000}` (divisores de 20000).
- **Fatorações**: `potrf` (Cholesky) + `getrf_nopiv` (LU) — as duas válidas no
  PaRSEC+CUDA. QR está fora (resultado numérico inválido no PaRSEC+GPU).
- **Escalonadores**: StarPU `{dmda, dmdas}`, PaRSEC `{lfq, gd}`.
- **Precisão**: FP32 e FP64.
- **THREADS** por nó (reserva 1 core p/ a GPU): poti=19, tupi=23.
- E1 **não coleta trace** — só `time,gflops` (o que ranqueia o tile).
- Binding de CPU **ligado** (Linux nativo); desligado automaticamente só sob WSL.

## Como rodar E1 no PCAD

```bash
# 1. (do PC local) sincronize o repo
bash scripts/pcad/sync_to_pcad.sh

# 2. (no front-end gppd-hpc) regenere os CSVs se preciso
nix develop --command Rscript final/doe_gpu_tile.r

# 3. submeta um job por partição
sbatch final/gpu_tile_poti.slurm
sbatch final/gpu_tile_tupi.slurm

# 4. (do PC local) puxe os resultados
bash scripts/pcad/get_from_pcad.sh
```

Resultados em `data/gpu_tile_{poti,tupi}_<jobid>/results.csv` (+ modelo StarPU
calibrado em `.starpu/`). Depois de E1, escolha o `b` campeão e gere E2 com
`final/doe_gpu_sched.r`.

## Smoke test local

```bash
DESIGN_FILE=final/doe_gpu_tile_poti.csv THREADS=8 \
  RESULTS_DIR=data/gpu_tile_local bash final/run_gpu_sched.sh
```
