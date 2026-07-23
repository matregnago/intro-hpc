library(tidyverse)

n <- 40000
b <- 1000

design <- tribble(
  ~precision, ~algorithm,     ~n, ~b, ~runtime, ~scheduler, ~rep,
  "FP64",     "potrf",        n,  b,  "starpu", "dmda",     1,
  "FP64",     "potrf",        n,  b,  "starpu", "dmdas",    1,
  "FP64",     "potrf",        n,  b,  "parsec", "gd",       1,
  "FP64",     "potrf",        n,  b,  "parsec", "lfq",      1,
  "FP64",     "getrf_nopiv",  n,  b,  "starpu", "dmda",     1,
  "FP64",     "getrf_nopiv",  n,  b,  "starpu", "dmdas",    1,
  "FP64",     "getrf_nopiv",  n,  b,  "parsec", "gd",       1,
  "FP64",     "getrf_nopiv",  n,  b,  "parsec", "lfq",      1,
)

write_csv(design, "doe_traces_tupi.csv", progress = FALSE)
