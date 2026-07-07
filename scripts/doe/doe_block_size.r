library(DoE.base)
library(tidyverse)

# E1 -- Selecao de tile (block size) na GPU. Objetivo: tracar a "montanha" de
# GFLOPS vs block size e LER O PICO = block size ideal, que ancora o `b` usado
# em E2b (doe_n_size_*.r) e E4 (doe_traces_*.r).
#
# N POR PRECISAO. FP32 N=80000 (divisores 3200/4000/5000 caem todos no topo da
# montanha); FP64 N=40000 (a FP64 da 4070 e capada, ~1 TFLOPS -- crescer N so
# multiplica o tempo sem revelar paralelismo novo). RAM nunca e gargalo (poti
# tem 96GB DDR5; a GPU so faz streaming de tiles).
#
# `b` DENSO onde o pico deve estar, e SEMPRE divisor de N (PaRSEC+CUDA da
# SIGSEGV/heap-corrupt em tile nao-divisor):
#   FP32: {1000,2000,3200,4000,5000,8000} (N/b 80->10)
#   FP64: {250,320,400,500,625,800,1000,2000,4000} (N/b 160->10, adensado
#     abaixo de 1000 -- e onde o pico FP64 cai)
#
# QR fica fora (dgeqrf/sgeqrf da resultado numerico invalido no PaRSEC+GPU).
# Escalonadores: StarPU {dmda,dmdas}, PaRSEC {lfq,gd}. N e tiles sao os MESMOS
# nos dois nos poti/tupi -- o que muda por no e so THREADS (env do
# slurm/block_size_{poti,tupi}.slurm), entao um unico doe_block_size.csv serve
# para os dois. Ordem das linhas randomizada (nao correlacionar deriva termica
# com nenhum fator).
#
# Custo estimado (reps=5): ~600 execucoes/no (FP32 240 + FP64 360). Tunavel por
# env: REPS (default 5).
#   Rscript scripts/doe/doe_block_size.r
# Saida: doe_block_size.csv

reps <- as.integer(Sys.getenv("REPS", "5"))

fator_algorithm     <- c("potrf", "getrf_nopiv")        # Cholesky + LU (QR fora)
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

precision_cfg <- list(
  FP32 = list(n = 80000L, blocks = c(1000, 2000, 3200, 4000, 5000, 8000)),
  FP64 = list(n = 40000L, blocks = c(250, 320, 400, 500, 625, 800, 1000, 2000, 4000))
)

make_design <- function(precision, n_fixed, blocks) {
  fac.design(
    nfactors = 3,
    replications = reps,
    repeat.only = FALSE,
    randomize = TRUE,
    seed = 1,
    nlevels = c(
      length(fator_algorithm),
      length(blocks),
      length(fator_runtime_sched)
    ),
    factor.names = list(
      algorithm = fator_algorithm,
      b = blocks,
      runtime_sched = fator_runtime_sched
    )
  ) |>
    as_tibble() |>
    mutate(
      algorithm = as.character(algorithm),
      b = as.integer(as.character(b)),
      runtime_sched = as.character(runtime_sched),
      precision = precision,
      n = n_fixed
    ) |>
    separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
    group_by(precision, algorithm, b, runtime, scheduler) |>
    mutate(rep = row_number()) |>
    ungroup() |>
    select(precision, algorithm, n, b, runtime, scheduler, rep)
}

design <- imap(precision_cfg, ~ make_design(.y, .x$n, .x$blocks)) |>
  bind_rows()

# Guarda dura: b TEM que dividir n (senao PaRSEC+CUDA quebra). Falha cedo.
stopifnot(all(design$n %% design$b == 0))

set.seed(1)
design <- slice_sample(design, prop = 1)

write_csv(design, "doe_block_size.csv", progress = FALSE)
cat(sprintf(
  "doe_block_size.csv : %d execucoes\n  FP32: N=%d, b={%s}\n  FP64: N=%d, b={%s}\n  reps=%d, schedulers={%s}, algos={%s}\n",
  nrow(design),
  precision_cfg$FP32$n, paste(precision_cfg$FP32$blocks, collapse = ","),
  precision_cfg$FP64$n, paste(precision_cfg$FP64$blocks, collapse = ","),
  reps, paste(fator_runtime_sched, collapse = ","),
  paste(fator_algorithm, collapse = ",")
))
