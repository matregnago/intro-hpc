library(tidyverse)

# E4-ideal -- coleta COMPLETA de rastro StarVZ (trace + DAG inline) de UMA
# execucao por celula: runtime x escalonador x algoritmo x precisao, no TAMANHO
# IDEAL de tile de cada (precisao,algoritmo) medido no E1b
# (job gpu_tile_ext_poti_797989: StarPU completo + PaRSEC nos divisores).
#
# Consumido por final/gpu_doe_sweep.sh via os SLURM gpu_trace_ideal_<node>.slurm
# com TRACE=1 TRACE_FULL=1 TRACE_DAG=1 TRACE_STATS=1 -> gera, por execucao:
#   StarPU: FxT (.fxt) + STARPU_BUS/WORKER_STATS no log;
#   PaRSEC: .prof (PINS task_profiler) + DAG .dot (CHAMELEON_PARSEC_DOT) +
#           device_show_statistics no log.
# Tudo o que o pipeline StarVZ precisa (paineis Espaco-Tempo; e p/ PaRSEC os
# converters parsec_dag_to_parquet.r / parsec_tasks_to_parquet.r -> dag/tasks/
# states.parquet + reconstrucao do Ready).
#
# Celulas = 4 (runtime,sched) x 2 algos x 2 precisoes = 16 por maquina, 1 rep
# (rastro e caro; 1 execucao basta p/ os paineis).
#
# N por precisao = o mesmo do tiling E1 (final/doe_gpu_tile.r): FP32 N=80000,
# FP64 N=40000, p/ capturar o rastro no regime exato em que E1 escolhe o tile.
# b SEMPRE divisor de N (PaRSEC+CUDA crasha em tile nao-divisor -- job 797989:
# 1500/3000/3500 deram 100% SIGSEGV/heap-corrupt no PaRSEC).
#
# Tile ideal por (precisao,algoritmo) = argmax GFlops sobre divisores do E1
# (StarPU); PaRSEC bate no FP32 e fica dentro do ruido no FP64. Os b abaixo sao
# os candidatos do E1b antigo -- CONFIRME com o pico do NOVO E1 (FP32 em N=80000)
# antes de submeter:
#   FP32 potrf       -> b=4000  (N=80000: 20x20 tiles, ~1.5k tarefas)
#   FP32 getrf_nopiv -> b=2500  (N=80000: 32x32, ~11k tarefas)
#   FP64 potrf       -> b=1000  (N=40000: 40x40, ~11k tarefas: FP64 e plato chato
#                       em [500,1000]; b=1000 e o pico do StarPU e mantem o rastro
#                       tratavel -- b menor explode p/ dezenas de milhar de tarefas)
#   FP64 getrf_nopiv -> b=1000  (idem)
# tupi reaproveita os tiles do poti: o run anterior da tupi (797982) teve a GPU
# disputada (StarPU OOM na calibracao + GFlops abaixo do 4070) -> sem dado
# confiavel de tile; reajustar quando houver um run limpo na tupi.
#
# Tunavel por env: REPS (replicacoes). N e b sao por precisao na tabela ideal_b.
#   Rscript final/doe_gpu_trace_ideal.r

reps <- as.integer(Sys.getenv("REPS", "1"))

# N por precisao espelha o tiling E1 (final/doe_gpu_tile.r): FP32 N=80000, FP64
# N=40000. b SEMPRE divisor do N da linha (guarda dura mais abaixo). Os b sao
# candidatos a pico do E1b antigo -- atualize com o pico do novo E1.
ideal_b <- tribble(
  ~precision, ~algorithm,        ~n,     ~b,
  "FP32",     "potrf",       80000L,  4000L,
  "FP32",     "getrf_nopiv", 80000L,  2500L,
  "FP64",     "potrf",       40000L,  1000L,
  "FP64",     "getrf_nopiv", 40000L,  1000L
)

runtime_sched <- tribble(
  ~runtime, ~scheduler,
  "starpu", "dmda",
  "starpu", "dmdas",
  "parsec", "lfq",
  "parsec", "gd"
)

design <- ideal_b |>
  crossing(runtime_sched) |>
  crossing(rep = seq_len(reps)) |>
  select(precision, algorithm, n, b, runtime, scheduler, rep) |>
  arrange(precision, algorithm, runtime, scheduler, rep)

# Guarda dura: b TEM que dividir n (senao PaRSEC+CUDA quebra). Falha cedo.
stopifnot(all(design$n %% design$b == 0))

for (node in c("poti", "tupi")) {
  out <- sprintf("final/doe_gpu_trace_ideal_%s.csv", node)
  write_csv(design, out, progress = FALSE)
  cat(sprintf("%-4s -> %s : %d execucoes (FP32 N=80000, FP64 N=40000, reps=%d)\n",
              node, out, nrow(design), reps))
}
