library(DoE.base)
library(tidyverse)

# E2 -- Impacto do escalonador na GPU: StarPU vs PaRSEC.
#
# Espelha a etapa "nmad_times" do Schnorr (WAMCA 2025): fixa N e o tile,
# cruza scheduler x runtime x algoritmo, randomizado e replicado. O estudo e
# somente FP64 (a RTX 4070 e hardware-capada em ~1 TFLOPS duplo, o que
# evidencia o custo do offload GPU vs CPU -- ponto central da comparacao).
#
# Tunavel por env: N (tamanho da matriz), B (tile), REPS (replicacoes).
#   Rscript scripts/doe/doe_gpu_sched.r

n_fixed <- as.integer(Sys.getenv("N", "19200"))
b_fixed <- as.integer(Sys.getenv("B", "960"))
reps <- as.integer(Sys.getenv("REPS", "10"))

fator_algorithm <- c("potrf", "geqrf")
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

design <- fac.design(
  nfactors = 2,
  replications = reps,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algorithm),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm = fator_algorithm,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm = as.character(algorithm),
    runtime_sched = as.character(runtime_sched)
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  mutate(precision = "FP64") |>
  group_by(precision, algorithm, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  mutate(n = n_fixed, b = b_fixed) |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

write_csv(design, "scripts/doe/doe_gpu_sched.csv", progress = FALSE)
print(design)
cat(sprintf("\n%d execucoes (N=%d, b=%d, reps=%d)\n",
            nrow(design), n_fixed, b_fixed, reps))
