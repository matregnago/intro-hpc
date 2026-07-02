library(DoE.base)
library(tidyverse)

# E1b -- selecao de tile (block size) na GPU, varredura ESTENDIDA em N maior.
#
# Motivacao: no E1 (final/doe_gpu_tile.r, N=20000) o Cholesky FP32 (spotrf) NAO
# atingiu o plato -- o pico de GFlops ficou na BORDA da grade de tiles (b=2000),
# ainda subindo. So estender o tile em N=20000 nao resolve: em N=20000 um
# b=4000 deixa a grade em 5x5=25 tiles, e a queda passa a confundir
# "tile grande demais p/ o kernel" com "ficou sem paralelismo". A solucao limpa
# e SUBIR O N: em N=40000 um b=4000 ainda da grade 10x10=100 tiles (mesma
# paralelizacao que o b=2000 tinha em N=20000), entao a curva vs b isola o
# efeito do tile. Aqui rodamos o FATORIAL COMPLETO (as 4 ops validas no
# PaRSEC+CUDA x 4 escalonadores) p/ ver como o tile otimo de CADA caso desloca
# com N -- nao so o spotrf.
#
# Diferencas vs doe_gpu_tile.r:
#   - N = 40000 (era 20000).
#   - b RESTRITO a DIVISORES de N=40000. O job 797989 mostrou que o caminho
#     Chameleon+PaRSEC+CUDA NAO trata o tile de resto: tiles nao-divisores
#     (1500/3000/3500) crasharam 100% no PaRSEC (SIGSEGV exit 139 / corrupcao de
#     heap glibc exit 134), enquanto o StarPU rodou. Para a comparacao ser justa
#     usamos o MESMO grid de tiles nos dois runtimes, e o unico grid seguro p/ o
#     PaRSEC e o de divisores. Em N=40000 (=2^6*5^4) os divisores em [400,4000]
#     sao {400,500,625,800,1000,1250,1600,2000,2500,4000}. ATENCAO: entre 2500 e
#     4000 nao existe divisor, entao a "montanha" do FP32 fica menos densa la em
#     cima do que no plano original (era esse o motivo do passo 500); aumentar a
#     densidade no topo exigiria um N muito maior (LCM dos tiles), fora do escopo.
#   - escreve final/doe_gpu_tile_ext_<node>.csv.
#
# Memoria (host-resident; a GPU so transita tiles): FP64 N=40000 = 12.8 GB,
# FP32 = 6.4 GB -- cabe folgado na RAM do host e na RTX 4070 (12 GB).
#
# Tunavel por env: REPS (replicacoes).
#   Rscript final/doe_gpu_tile_ext.r

reps <- as.integer(Sys.getenv("REPS", "3"))

fator_precision <- c("FP32", "FP64")
fator_algorithm <- c("potrf", "getrf_nopiv")           # Cholesky + LU (QR fora)
fator_block <- c(400, 500, 625, 800, 1000,
                 1250, 1600, 2000, 2500, 4000)         # so divisores de N=40000
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

# Config por no. N igual nos dois (comparacao 4070 vs 4090 no mesmo problema);
# `threads` e so documentacao (vai como THREADS no runner, nao como coluna).
nodes <- list(
  poti = list(n = 40000L, threads = 19L),
  tupi = list(n = 40000L, threads = 23L)
)

make_design <- function(n_fixed) {
  fac.design(
    nfactors = 4,
    replications = reps,
    repeat.only = FALSE,
    randomize = TRUE,
    seed = 1,
    nlevels = c(
      length(fator_precision),
      length(fator_algorithm),
      length(fator_block),
      length(fator_runtime_sched)
    ),
    factor.names = list(
      precision = fator_precision,
      algorithm = fator_algorithm,
      b = fator_block,
      runtime_sched = fator_runtime_sched
    )
  ) |>
    as_tibble() |>
    mutate(
      precision = as.character(precision),
      algorithm = as.character(algorithm),
      b = as.integer(as.character(b)),
      runtime_sched = as.character(runtime_sched)
    ) |>
    separate(runtime_sched, into = c("runtime", "scheduler"), sep = ":") |>
    group_by(precision, algorithm, b, runtime, scheduler) |>
    mutate(rep = row_number()) |>
    ungroup() |>
    mutate(n = n_fixed) |>
    select(precision, algorithm, n, b, runtime, scheduler, rep)
}

for (node in names(nodes)) {
  cfg <- nodes[[node]]
  design <- make_design(cfg$n)
  out <- sprintf("final/doe_gpu_tile_ext_%s.csv", node)
  write_csv(design, out, progress = FALSE)
  cat(sprintf(
    "%-4s -> %s : %d execucoes (N=%d, b={%s}, reps=%d, THREADS=%d)\n",
    node, out, nrow(design), cfg$n,
    paste(fator_block, collapse = ","), reps, cfg$threads
  ))
}
