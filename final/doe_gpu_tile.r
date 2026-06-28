library(DoE.base)
library(tidyverse)

# E1 -- Selecao de tile (block size) na GPU, gerada POR NO (poti e tupi).
#
# Espelha a etapa "block-size" do Schnorr (WAMCA 2025): fixa N, varia o tile `b`
# e cruza com (precisao x algoritmo x escalonador). A saida ancora o `b` usado
# em E2 (doe_gpu_sched.r). Diferente da versao antiga, varre os QUATRO
# escalonadores do estudo (nao so os "campeoes") e usa as DUAS fatoracoes
# validas no PaRSEC+CUDA -- Cholesky (potrf) e LU (getrf_nopiv). QR fica de fora:
# dgeqrf/sgeqrf dao resultado numerico INVALIDO no PaRSEC+GPU.
#
# Dois arquivos sao escritos, um por no:
#   scripts/doe/doe_gpu_tile_poti.csv   (i7-14700KF 20 cores, RTX 4070 12GB)
#   scripts/doe/doe_gpu_tile_tupi.csv   (i9-14900KF 24 cores, RTX 4090 24GB)
# N e os tiles sao IGUAIS nos dois (comparacao 4070 vs 4090 no mesmo problema);
# o que muda por no e o THREADS passado ao runner (reserva 1 core p/ dirigir a
# GPU): poti=19, tupi=23. Esse valor NAO vai no CSV (gpu_doe_sweep.sh le THREADS
# do ambiente) -- esta documentado aqui e impresso ao final.
#
# Tunavel por env: REPS (replicacoes).
#   Rscript scripts/doe/doe_gpu_tile.r

reps <- as.integer(Sys.getenv("REPS", "3"))

fator_precision <- c("FP32", "FP64")
fator_algorithm <- c("potrf", "getrf_nopiv")           # Cholesky + LU (QR fora)
fator_block <- c(400, 500, 800, 1000, 2000)            # divisores de 20000
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

# Config por no. N igual nos dois; `threads` e so documentacao (vai como
# THREADS no runner, nao como coluna do CSV).
nodes <- list(
  poti = list(n = 20000L, threads = 19L),
  tupi = list(n = 20000L, threads = 23L)
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
  out <- sprintf("final/doe_gpu_tile_%s.csv", node)
  write_csv(design, out, progress = FALSE)
  cat(sprintf(
    "%-4s -> %s : %d execucoes (N=%d, b={%s}, reps=%d, THREADS=%d)\n",
    node, out, nrow(design), cfg$n,
    paste(fator_block, collapse = ","), reps, cfg$threads
  ))
}
