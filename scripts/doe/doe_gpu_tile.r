library(DoE.base)
library(tidyverse)

# E1 -- Selecao de tile (block size) na GPU.
#
# Espelha a etapa "block-size" do Schnorr (WAMCA 2025): fixa N e o escalonador
# campeao de cada runtime (StarPU=dmda, PaRSEC=lfq), varia apenas o tile `b`,
# por (runtime x algoritmo x precisao). A saida ancora o `b` usado em E2.
#
# Tunavel por env: N (tamanho da matriz), REPS (replicacoes).
#   Rscript scripts/doe/doe_gpu_tile.r

n_fixed <- as.integer(Sys.getenv("N", "19200"))
reps <- as.integer(Sys.getenv("REPS", "5"))

fator_precision <- c("FP32", "FP64")
fator_algorithm <- c("potrf", "geqrf")
fator_block <- c(480, 720, 960, 1440, 1920)
fator_runtime_sched <- c("starpu:dmda", "parsec:lfq")

design <- fac.design(
  nfactors = 4,
  replications = reps,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_precision),
    length(fator_algorithm),
    length(fator_block),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    precision = fator_precision,
    algorithm = fator_algorithm,
    b = fator_block,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    precision = as.character(precision),
    algorithm = as.character(algorithm),
    b = as.integer(as.character(b)),
    runtime_sched = as.character(runtime_sched)
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  group_by(precision, algorithm, b, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  mutate(n = n_fixed) |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

write_csv(design, "scripts/doe/doe_gpu_tile.csv", progress = FALSE)
print(design)
cat(sprintf("\n%d execucoes (N=%d, reps=%d)\n", nrow(design), n_fixed, reps))
