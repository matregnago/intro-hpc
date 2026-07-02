library(DoE.base)
library(tidyverse)

# E2b -- Escalonamento vs tamanho do problema (N) na GPU, no pico de tile (E1).
# Objetivo: tracar GFLOPS vs N e ler a ASSIMPTOTA de saturacao da RTX 4070,
# ancorada no block size IDEAL (b=500) achado em E1 por algoritmo/escalonador
# (doe_gpu_tile.r -> gpu_tile_poti_798331: pico FP64 quase sempre em b=500).
#
# FP64-ONLY (a motivacao do trabalho e densa em precisao dupla);
#   o pico de tile FP64 foi determinado na poti, entao so poti aqui.
# b=500 fixo e sempre divisor de N (PaRSEC+CUDA quebra com SIGSEGV se b nao
#   divide n -- job 797989). Traca QR de novo fora (dgeqrf invalido no PaRSEC+GPU).
#
# N varia em multiplos de 500, mantendo N/b entre 10 e 80 (abaixo de 10 a GPU
# nao satura e o erro numerico cresce; acima de 40000 a memoria de 12GB doi):
#   {5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000}
# Sobrescreva via env: N_SIZES="5000 10000 20000 40000" REPS=3 Rscript ... .
#
# Custo (reps=5, 8 Ns): 8 N x 2 alg x 4 sched = 320 execucoes (160/starpu +
# 160/parsec), + 2 passes de calibracao StarPU por (alg,N) na fase starpu.
#   Rscript final/doe_gpu_n_size_poti.r
#
# Escreve final/doe_gpu_n_size_poti.csv.

reps <- as.integer(Sys.getenv("REPS", "5"))
n_sizes <- as.integer(strsplit(
  trimws(Sys.getenv("N_SIZES", "5000 10000 15000 20000 25000 30000 35000 40000")),
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

# Guarda dura: b tem que dividir todos os N (senao PaRSEC+CUDA quebra).
stopifnot(all(design$n %% design$b == 0))

# Randomiza a ordem de execucao para nao correlacionar deriva termica com N.
set.seed(1)
design <- slice_sample(design, prop = 1)

out <- "final/doe_gpu_n_size_poti.csv"
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
