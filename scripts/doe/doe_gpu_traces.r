library(DoE.base)
library(tidyverse)

# E4 -- Rastros (StarVZ) para explicar o *porque* do ranking de escalonadores.
#
# Subconjunto representativo da E2 (doe_gpu_sched.r) com o MESMO formato de CSV
# que o scripts/gpu_doe_sweep.sh consome, mas rodado com TRACE=1 pelo SLURM para
# capturar os arquivos de rastro (StarPU FxT / PaRSEC profile) por execucao.
# Cruza scheduler x runtime x algoritmo x precisao nas configs representativas
# (melhor vs pior escalonador de cada runtime), com poucas replicacoes -- rastro
# e caro e basta 1-2 execucoes por celula para os paineis Espaco-Tempo.
#
# N e b ancorados em E1/E2 (b favorece os tiles do cuBLAS na GPU).
# Tunavel por env: N (matriz), B (tile), REPS (replicacoes).
#   Rscript scripts/doe/doe_gpu_traces.r   # -> scripts/doe/doe_gpu_traces.csv

n_fixed <- as.integer(Sys.getenv("N", "19200"))
b_fixed <- as.integer(Sys.getenv("B", "960"))
reps <- as.integer(Sys.getenv("REPS", "1"))

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

write_csv(design, "scripts/doe/doe_gpu_traces.csv", progress = FALSE)
print(design)
cat(sprintf("\n%d execucoes com rastro (N=%d, b=%d, reps=%d)\n",
            nrow(design), n_fixed, b_fixed, reps))
