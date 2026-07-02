library(DoE.base)
library(tidyverse)

# E2b (tupi) -- Escalonamento vs tamanho do problema (N) na GPU, no pico de
# tile (E1) da tupi. Gegenparte de final/doe_gpu_n_size_poti.r mas com o pico
# de tile da RTX 4090: b=1000 (do job gpu_tile_tupi_20260701_095748: StarPU
# pica em b=1000, PaRSEC em b=2000 -- usamos 1000 como compromisso;
# 1000 tambem divide todos os N do sweep).
#
# FP64-ONLY; mesmos Ns e estrutura do DoE poti, so troca b=500 -> b=1000.
# b=1000 fixo e sempre divisor de N (PaRSEC+CUDA quebra com SIGSEGV se b nao
#   divide n). Traca QR fora (dgeqrf invalido no PaRSEC+GPU).
#
# N varia em multiplos de 1000, mantendo N/b >= 5 (abaixo disso a GPU nao
# satura e o erro numerico cresce). ESTENDIDO ate 60000 depois que o job
# gpu_n_size_tupi_20260701_115305 mostrou a curva ainda subindo em N=40000
# (dgetrf 1517 -> 1582 no ultimo degrau, sem plato):
#   {5000..40000 passo 5000} + {50000, 60000}
# N=50000 (20 GB) ja NAO cabe nos 21.2 GiB da RTX 4090 -> o regime out-of-core
# (evicao, como no poti FP64) entra de proposito no fim da curva.
# Sobrescreva via env: N_SIZES="5000 10000 20000 40000" REPS=5 Rscript ... .
#
# Custo (reps=3, 10 Ns, tempos medidos do job 20260701_115305 + escala N^3):
# 10 N x 2 alg x 4 sched x 3 reps = 240 execucoes ~= 67 min + ~11 min de
# calibracao StarPU (2 passes por (alg,N)) + overhead de shell/rsync ->
# ~1h40; alocar >= 3h de folga.
#   Rscript final/doe_gpu_n_size_tupi.r
#
# Escreve final/doe_gpu_n_size_tupi.csv.

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

# Guarda dura: b tem que dividir todos os N (senao PaRSEC+CUDA quebra).
stopifnot(all(design$n %% design$b == 0))

# Randomiza a ordem de execucao para nao correlacionar deriva termica com N.
set.seed(1)
design <- slice_sample(design, prop = 1)

out <- Sys.getenv("OUT", "final/doe_gpu_n_size_tupi.csv")
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
