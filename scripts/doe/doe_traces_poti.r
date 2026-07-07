library(tidyverse)

# E4 (poti) -- DoE da captura de rastros (StarVZ) na RTX 4070: b=500 (pico de
# E1 em FP64), N=40000, 1 rep -- o objetivo aqui e o TRACE em si (space-time,
# k-iteration, task times, transferencias), nao a estatistica de tempo, entao
# nao ha replicacao. Fatoracoes: potrf + getrf_nopiv (QR fora). Escalonadores:
# StarPU {dmda,dmdas}, PaRSEC {gd,lfq} -> 2 x 4 x 1 = 8 execucoes.
#
# Roda via slurm/traces_poti.slurm, que exporta TRACE=1 TRACE_FULL=1
# TRACE_DAG=0 TRACE_STATS=1 (dag fica pra um job dedicado -- nesses N o
# grapher polui o tempo medido) e CALIB_PASSES=2 (calibra dmda/dmdas antes das
# medidas cronometradas, scripts/run.sh::calibrate_starpu).
#
#   Rscript scripts/doe/doe_traces_poti.r
# Saida: doe_traces_poti.csv (8 linhas)

n <- 40000L
b <- 500L
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

write_csv(design, "doe_traces_poti.csv", progress = FALSE)
message("escrito doe_traces_poti.csv (", nrow(design), " linhas)")
