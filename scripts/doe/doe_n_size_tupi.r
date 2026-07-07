library(DoE.base)
library(tidyverse)

# E2b (tupi) -- contraparte de doe_n_size_poti.r no pico de tile da RTX 4090
# (b=1000 em FP64 -- ver doe_block_size.r / gflops_vs_b_compare.png).
#
# FP64-ONLY; mesma estrutura do DoE poti, so troca b=500 -> b=1000. N varia em
# multiplos de 1000, mantendo N/b >= 5 (abaixo disso a GPU nao satura e o erro
# numerico cresce). ESTENDIDO ate 60000 (a curva ainda subia em N=40000 sem
# platô):
#   {5000..40000 passo 5000} + {50000, 60000}
# N=50000 (20 GB) ja nao cabe nos 21.2 GiB uteis da RTX 4090 -> o regime
# out-of-core (evicao, como no poti FP64) entra de proposito no fim da curva.
# Sobrescreva via env: N_SIZES="5000 10000 20000 40000" REPS=5 Rscript ... .
#
# Fatoracoes: potrf + getrf_nopiv (QR fora). Escalonadores: StarPU
# {dmda,dmdas}, PaRSEC {lfq,gd}.
#
# Custo (reps=3, 10 Ns): 10 N x 2 alg x 4 sched x 3 reps = 240 execucoes.
#   Rscript scripts/doe/doe_n_size_tupi.r
# Saida: doe_n_size_tupi.csv

reps <- as.integer(Sys.getenv("REPS", "3"))
n_sizes <- as.integer(strsplit(
  trimws(Sys.getenv(
    "N_SIZES",
    "5000 10000 15000 20000 25000 30000 35000 40000 50000 60000")),
  "\\s+")[[1]])
block <- as.integer(Sys.getenv("BLOCK", "1000"))

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

out <- "doe_n_size_tupi.csv"
write_csv(design, out, progress = FALSE)

cat(sprintf(
  "tupi -> %s : %d execucoes (THREADS=23, b=%d, reps=%d)\n",
  out, nrow(design), block, reps
))
cat(sprintf(
  "  N={%s}\n  schedulers={%s}\n  algos={%s}\n",
  paste(n_sizes, collapse = ","), paste(fator_runtime_sched, collapse = ","),
  paste(fator_algorithm, collapse = ",")
))
