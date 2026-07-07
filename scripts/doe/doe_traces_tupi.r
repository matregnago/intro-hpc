library(tidyverse)

# E4 (tupi) -- contraparte de doe_traces_poti.r na RTX 4090: b=1000 (pico de
# E1 em FP64), N=40000, 1 rep. Fatoracoes: potrf + getrf_nopiv (QR fora).
# Escalonadores: StarPU {dmda,dmdas}, PaRSEC {gd,lfq} -> 2 x 4 x 1 = 8 execucoes.
#
# Roda via slurm/traces_tupi.slurm (TRACE=1 TRACE_FULL=1 TRACE_DAG=0
# TRACE_STATS=1, CALIB_PASSES=2 -- ver doe_traces_poti.r).
#
#   Rscript scripts/doe/doe_traces_tupi.r
# Saida: doe_traces_tupi.csv (8 linhas)

n <- 40000L
b <- 1000L
stopifnot(n %% b == 0)

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
message("escrito doe_traces_tupi.csv (", nrow(design), " linhas)")
