library(DoE.base)
library(tidyverse)

# E1 -- Selecao de tile (block size) na GPU. Objetivo: tracar a "montanha" de
# GFLOPS vs block size e LER O PICO = block size ideal, que ancora o `b` usado
# em E2 (doe_gpu_sched.r).
#
# Reescrito para problemas que SATURAM a GPU (RTX 4070), depois de constatar nos
# dados que N pequeno deixava as execucoes curtas demais (FP32 spotrf em N=40000
# dava ~1.6s) e a curva GFLOPS-vs-b nunca chegava no pico. SUPERSEDE a antiga
# varredura ESTENDIDA (doe_gpu_tile_ext.r, N=40000): la os divisores de 40000
# deixavam um BURACO entre 2500 e 4000, ralando o topo da montanha FP32. Em
# N=80000 (=2^7*5^4) os divisores 3200/4000/5000 caem todos no topo, densos.
# Mudancas:
#
#   * N POR PRECISAO. FP32 vai a N=80000 (steady-state ~11-22s StarPU, ~33-66s
#     PaRSEC); FP64 fica em N=40000 -- la ja dura 20-47s e a FP64 da 4070 e
#     capada (~1 TFLOPS), entao crescer N so multiplica o tempo (cresce com N^3)
#     sem revelar paralelismo novo. RAM nunca e gargalo: poti tem 96GB DDR5
#     (FP32 N=80000 = 25.6GB no host; a GPU so faz streaming de tiles).
#
#   * b DENSO onde o pico deve estar, e SEMPRE divisor de N (PaRSEC+CUDA quebra
#     com SIGSEGV se b nao divide n -- job 797989: 1500/3000/3500 deram 100%
#     SIGSEGV/heap-corrupt no PaRSEC). FP32: {1000,2000,3200,4000,5000,8000}
#     (N/b de 80 a 10) -- adensado em 3200-5000, a regiao provavel do pico, com
#     as asas em 1000/2000 (subida) e 8000 (descida, N/b=10). FP64:
#     {250,320,400,500,625,800,1000,2000,4000} (N/b de 160 a 10) -- adensado
#     ABAIXO de 1000 porque a varredura antiga {1000,2000,4000} achava o pico
#     SEMPRE em b=1000, a borda esquerda da grade: a curva FP64 NAO
#     e plana, ela sobe ate b pequeno e o topo da montanha ficava fora da grade.
#     Os 6 blocos novos sao todos divisores de 40000 (=2^6*5^4) abaixo de 1000.
#
#   * QR continua FORA: dgeqrf/sgeqrf dao resultado numerico INVALIDO no
#     PaRSEC+GPU. Varre os QUATRO escalonadores (StarPU dmda/dmdas, PaRSEC
#     lfq/gd) -- nenhum e GPU-aware no sentido do dmda; a escolha CPU-vs-GPU do
#     PaRSEC esta no device layer (parsec_gpu_get_best_device), igual p/ lfq/gd.
#
# Dois arquivos sao escritos, um por no (conteudo identico; N e tiles iguais nos
# dois para comparar 4070 vs 4090 no mesmo problema; o que muda por no e o
# THREADS que o runner le do ambiente, documentado abaixo e impresso ao final):
#   final/doe_gpu_tile_poti.csv   (i7-14700KF 20c, RTX 4070 12GB, THREADS=19)
#   final/doe_gpu_tile_tupi.csv   (i9-14900KF 24c, RTX 4090 24GB, THREADS=23)
#
# Sem traces nesta fase (so time/gflops). A ordem das linhas e randomizada para
# nao correlacionar deriva termica com nenhum fator.
#
# Custo estimado (reps=5): ~600 execucoes/no (FP32 240 + FP64 360). Tunavel por
# env: REPS (default 5).
#   Rscript final/doe_gpu_tile.r

reps <- as.integer(Sys.getenv("REPS", "5"))

fator_algorithm     <- c("potrf", "getrf_nopiv")        # Cholesky + LU (QR fora)
fator_runtime_sched <- c(
  "starpu:dmda",
  "starpu:dmdas",
  "parsec:lfq",
  "parsec:gd"
)

# N e blocos POR PRECISAO -- b SEMPRE divisor de N.
precision_cfg <- list(
  FP32 = list(n = 80000L, blocks = c(1000, 2000, 3200, 4000, 5000, 8000)),
  FP64 = list(n = 40000L, blocks = c(250, 320, 400, 500, 625, 800, 1000, 2000, 4000))
)

# THREADS por no (reserva 1 core p/ dirigir a GPU) -- documentacao; vai como env
# THREADS ao runner, nao como coluna do CSV.
nodes <- list(
  poti = list(threads = 19L),
  tupi = list(threads = 23L)
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

# Randomiza a ordem de execucao misturando precisao/algoritmo/tile.
set.seed(1)
design <- slice_sample(design, prop = 1)

for (node in names(nodes)) {
  cfg <- nodes[[node]]
  out <- sprintf("final/doe_gpu_tile_%s.csv", node)
  write_csv(design, out, progress = FALSE)
  cat(sprintf(
    "%-4s -> %s : %d execucoes (THREADS=%d)\n",
    node, out, nrow(design), cfg$threads
  ))
}

cat(sprintf(
  "  FP32: N=%d, b={%s}\n  FP64: N=%d, b={%s}\n  reps=%d, schedulers={%s}, algos={%s}\n",
  precision_cfg$FP32$n, paste(precision_cfg$FP32$blocks, collapse = ","),
  precision_cfg$FP64$n, paste(precision_cfg$FP64$blocks, collapse = ","),
  reps, paste(fator_runtime_sched, collapse = ","),
  paste(fator_algorithm, collapse = ",")
))
