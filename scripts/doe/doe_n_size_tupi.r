library(DoE.base)
library(tidyverse)

reps <- 3
n_sizes <- c(5000, 10000, 15000, 20000, 25000, 30000,
             35000, 40000, 50000, 60000)
block <- 1000

fator_algorithm     <- c("potrf", "getrf_nopiv")
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

set.seed(1)
design <- slice_sample(design, prop = 1)

out <- "doe_n_size_tupi.csv"
write_csv(design, out, progress = FALSE)
