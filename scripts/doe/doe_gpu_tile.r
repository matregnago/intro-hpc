library(DoE.base)
library(tidyverse)

# E1 -- Selecao de tile (block size) na GPU. Objetivo: tracar a "montanha" de
# GFLOPS vs block size e LER O PICO = block size ideal, que ancora o `b` usado
# em E2 (doe_gpu_sched.r).
#
# Reescrito para um problema que SATURA a GPU (RTX 4070), depois de constatar
# nos dados que N pequeno deixava as execucoes curtas demais e a curva
# GFLOPS-vs-b nunca chegava no pico (subia monotonicamente ate o maior b
# testado). O estudo e SOMENTE FP64 (a 4070 e capada em ~1 TFLOPS duplo, o que
# evidencia o custo do offload GPU vs CPU -- ponto central da comparacao).
#
#   * N=40000 -- la ja dura 20-47s e a FP64 da 4070 e capada (~1 TFLOPS), entao
#     crescer N so multiplica o tempo (cresce com N^3) sem revelar paralelismo
#     novo. RAM nunca e gargalo: poti tem 96GB DDR5 (FP64 N=40000 = 12.8GB no
#     host; a GPU so faz streaming de tiles).
#
#   * b DENSO onde o pico deve estar, e SEMPRE divisor de N (PaRSEC+CUDA quebra
#     com SIGSEGV se b nao divide n). {1000,2000,4000} (N/b de 40 a 10) -- a
#     curva FP64 e plana e as execucoes sao caras; b=8000 em N=40000 daria N/b=5
#     (paralelismo faminto, ponto nao-informativo), entao fica de fora.
#
#   * QR continua FORA: dgeqrf dao resultado numerico INVALIDO no PaRSEC+GPU.
#     Varre os QUATRO escalonadores do estudo (StarPU dmda/dmdas, PaRSEC lfq/gd)
#     -- nenhum e GPU-aware no sentido do dmda; a escolha CPU-vs-GPU do PaRSEC
#     esta no device layer (parsec_gpu_get_best_device), igual p/ lfq/gd.
#
# Sem traces nesta fase (so time/gflops). A ordem das linhas e randomizada para
# nao correlacionar deriva termica com nenhum fator.
#
# Custo estimado (reps=3): ~72 execucoes/no, ~45min de walltime -- pedir >=1h no
# SLURM. Tunavel por env: REPS (default 3).
#   Rscript scripts/doe/doe_gpu_tile.r

reps <- as.integer(Sys.getenv("REPS", "3"))

fator_algorithm     <- c("potrf", "getrf_nopiv")        # Cholesky + LU (QR fora)
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

# N e blocos -- b SEMPRE divisor de N.
n_fixed <- 40000L
blocks  <- c(1000, 2000, 4000)

# Config por no. N e tiles IGUAIS nos dois (comparacao 4070 vs 4090 no mesmo
# problema); o que muda e o THREADS passado ao runner (reserva 1 core p/ dirigir
# a GPU): poti=19, tupi=23. Esse valor NAO vai no CSV (gpu_doe_sweep.sh le
# THREADS do ambiente) -- esta documentado aqui e impresso ao final.
nodes <- list(
  poti = list(threads = 19L),
  tupi = list(threads = 23L)
)

design <- fac.design(
  nfactors = 3,
  replications = reps,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algorithm),
    length(blocks),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm = fator_algorithm,
    b = blocks,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm = as.character(algorithm),
    b = as.integer(as.character(b)),
    runtime_sched = as.character(runtime_sched),
    precision = "FP64",
    n = n_fixed
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  group_by(precision, algorithm, b, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

# Guarda dura: b TEM que dividir n (senao PaRSEC+CUDA quebra). Falha cedo.
stopifnot(all(design$n %% design$b == 0))

# Randomiza a ordem de execucao misturando algoritmo/tile.
set.seed(1)
design <- slice_sample(design, prop = 1)

for (node in names(nodes)) {
  cfg <- nodes[[node]]
  out <- sprintf("scripts/doe/doe_gpu_tile_%s.csv", node)
  write_csv(design, out, progress = FALSE)
  cat(sprintf(
    "%-4s -> %s : %d execucoes (THREADS=%d)\n",
    node, out, nrow(design), cfg$threads
  ))
}

cat(sprintf(
  "  FP64: N=%d, b={%s}\n  reps=%d, schedulers={%s}, algos={%s}\n",
  n_fixed, paste(blocks, collapse = ","),
  reps, paste(fator_runtime_sched, collapse = ","),
  paste(fator_algorithm, collapse = ",")
))
