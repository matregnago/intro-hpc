library(DoE.base)
library(tidyverse)

fator_n_size <- 19200
fator_algoritmos <- c("potrf", "geqrf")
fator_iterations <- c(96, 75, 60, 50, 40, 30, 20)
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
    length(fator_iterations),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm     = fator_algoritmos,
    iterations    = fator_iterations,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm     = as.character(algorithm),
    runtime_sched = as.character(runtime_sched),
    iterations    = as.integer(as.character(iterations))
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  mutate(
    n = fator_n_size,
    b = fator_n_size %/% iterations
  ) |>
  arrange(runtime) |>
  relocate(algorithm, n, b, runtime, scheduler, iterations) |>
  print()

write_csv(design, "doe_block_size.csv", progress = FALSE)
