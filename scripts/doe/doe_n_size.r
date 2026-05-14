library(DoE.base)
library(tidyverse)

bloco <- 480
fator_n_size <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
fator_algoritmos <- c("dpotrf", "dgeqrf")
fator_runtime_sched <- c(
  "starpu:lws",
  "starpu:ws",
  "parsec:lfq",
  "parsec:gd"
)

design <- fac.design(
  nfactors = 3,
  replications = 5,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algoritmos),
    length(fator_n_size),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm = fator_algoritmos,
    n = fator_n_size,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm = as.character(algorithm),
    runtime_sched = as.character(runtime_sched),
    n = as.integer(as.character(n))
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  mutate(
    b = bloco,
    n = n * bloco
  ) |>
  relocate(algorithm, n, b, runtime, scheduler) |>
  print()

write_csv(design, "doe_n_size.csv", progress = FALSE)
