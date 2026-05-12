library(DoE.base)
library(tidyverse)

n_size <- 19200
iterations <- 40
fator_algoritmos <- c("potrf", "geqrf")
fator_runtime_sched <- c(
  "starpu:lws",
  "starpu:ws",
  "parsec:lfq",
  "parsec:gd"
)

design <- fac.design(
  nfactors = 2,
  replications = 1,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algoritmos),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm     = fator_algoritmos,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm     = as.character(algorithm),
    runtime_sched = as.character(runtime_sched)
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  mutate(
    n = n_size,
    b = n_size %/% iterations
  ) |>
  arrange(runtime) |>
  relocate(algorithm, n, b, runtime, scheduler) |>
  print()

write_csv(design, "doe_traces.csv", progress = FALSE)
