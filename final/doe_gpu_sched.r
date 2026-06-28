library(DoE.base)
library(tidyverse)

# E2 -- Impacto do escalonador na GPU: StarPU vs PaRSEC.
#
# Espelha a etapa "nmad_times" do Schnorr (WAMCA 2025): fixa N e o tile,
# cruza scheduler x runtime x algoritmo x precisao, randomizado e replicado.
# Aqui o segundo fator (analogo a distribuicao PxQ do Schnorr, que so existe
# com MPI) e a PRECISAO {FP32, FP64} -- onde a GPU de consumo ganha (FP32) ou
# empata com a CPU (FP64).
#
# Tunavel por env: N (tamanho da matriz), B (tile), REPS (replicacoes).
#   Rscript scripts/doe/doe_gpu_sched.r

n_fixed <- as.integer(Sys.getenv("N", "19200"))
b_fixed <- as.integer(Sys.getenv("B", "960"))
reps <- as.integer(Sys.getenv("REPS", "10"))

fator_precision <- c("FP32", "FP64")
fator_algorithm <- c("potrf", "geqrf")
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
    length(fator_precision),
    length(fator_algorithm),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    precision = fator_precision,
    algorithm = fator_algorithm,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    precision = as.character(precision),
    algorithm = as.character(algorithm),
    runtime_sched = as.character(runtime_sched)
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  group_by(precision, algorithm, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  mutate(n = n_fixed, b = b_fixed) |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

write_csv(design, "final/doe_gpu_sched.csv", progress = FALSE)
print(design)
cat(sprintf("\n%d execucoes (N=%d, b=%d, reps=%d)\n",
            nrow(design), n_fixed, b_fixed, reps))
