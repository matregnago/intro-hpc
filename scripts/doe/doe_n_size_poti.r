library(DoE.base)
library(tidyverse)

# E2b (poti) -- Escalonamento vs tamanho do problema (N) na GPU, no pico de
# tile de E1 (b=500 em FP64 -- ver doe_block_size.r / gflops_vs_b_compare.png).
# Objetivo: tracar GFLOPS vs N e ler a ASSIMPTOTA de saturacao da RTX 4070.
#
# FP64-ONLY (foco do estudo). b=500 fixo, sempre divisor de N (PaRSEC+CUDA
# quebra com SIGSEGV se b nao divide n). N/b entre 10 e 120 (abaixo de 10 a
# GPU nao satura e o erro numerico cresce); ESTENDIDO ate 60000 (como na tupi
# -- ver doe_n_size_tupi.r) depois que a curva ainda subia em N=40000 sem
# platô:
#   {5000, 10000, 20000, 40000, 50000, 60000}
# N=50000/60000 (20/28.8 GB) ja nao cabem nos 12GB da RTX 4070 -> o regime
# out-of-core (evicao) entra de proposito no fim da curva.
# Sobrescreva via env: N_SIZES="5000 10000 20000 40000" REPS=3 Rscript ... .
#
# Fatoracoes: potrf + getrf_nopiv (QR fora -- resultado numerico invalido no
# PaRSEC+GPU). Escalonadores: StarPU {dmda,dmdas}, PaRSEC {lfq,gd}.
#
# Custo (reps=5, 6 Ns): 6 N x 2 alg x 4 sched x 5 reps = 240 execucoes
# (120/starpu + 120/parsec), + calibracao StarPU por (alg,N) na fase starpu.
#   Rscript scripts/doe/doe_n_size_poti.r
# Saida: doe_n_size_poti.csv

reps <- as.integer(Sys.getenv("REPS", "5"))
n_sizes <- as.integer(strsplit(
  trimws(Sys.getenv("N_SIZES", "5000 10000 20000 40000 50000 60000")),
  "\\s+")[[1]])
block <- as.integer(Sys.getenv("BLOCK", "500"))

fator_algorithm     <- c("potrf", "getrf_nopiv")        # Cholesky + LU (QR fora)
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

design <- fac.design(
  nfactors = 3,
  replications = reps,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algorithm),
    length(n_sizes),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm = fator_algorithm,
    n = n_sizes,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm = as.character(algorithm),
    n = as.integer(as.character(n)),
    runtime_sched = as.character(runtime_sched),
    precision = "FP64",
    b = block
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  group_by(precision, algorithm, n, b, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

stopifnot(all(design$n %% design$b == 0))

set.seed(1)
design <- slice_sample(design, prop = 1)

out <- "doe_n_size_poti.csv"
write_csv(design, out, progress = FALSE)

cat(sprintf(
  "poti -> %s : %d execucoes (THREADS=19, b=%d, reps=%d)\n",
  out, nrow(design), block, reps
))
cat(sprintf(
  "  N={%s}\n  schedulers={%s}\n  algos={%s}\n",
  paste(n_sizes, collapse = ","), paste(fator_runtime_sched, collapse = ","),
  paste(fator_algorithm, collapse = ",")
))
