library(DoE.base)
library(tidyverse)


reps <- 5

fator_algorithm     <- c("potrf", "getrf_nopiv")
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

n_fixed <- 40000
fator_b <- c(250, 320, 400, 500, 625, 800, 1000, 2000, 4000)

design <- fac.design(
  nfactors = 3,
  replications = reps,
  repeat.only = FALSE,
  randomize = TRUE,
  seed = 1,
  nlevels = c(
    length(fator_algorithm),
    length(fator_b),
    length(fator_runtime_sched)
  ),
  factor.names = list(
    algorithm = fator_algorithm,
    b = fator_b,
    runtime_sched = fator_runtime_sched
  )
) |>
  as_tibble() |>
  mutate(
    algorithm = as.character(algorithm),
    b = as.integer(as.character(b)),
    runtime_sched = as.character(runtime_sched),
    precision = "FP64",
    n = n_fixed
  ) |>
  separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
  group_by(precision, algorithm, b, runtime, scheduler) |>
  mutate(rep = row_number()) |>
  ungroup() |>
  select(precision, algorithm, n, b, runtime, scheduler, rep)

set.seed(1)
design <- slice_sample(design, prop = 1)

write_csv(design, "doe_block_size.csv", progress = FALSE)
