#!/usr/bin/env Rscript

# E1 -- "Montanha" de GFLOPS vs block size, para LER O PICO = block size ideal
# (que ancora o `b` de E2, doe_gpu_sched.r). Le o sweep de tile produzido por
# scripts/doe/doe_gpu_tile.r -> gpu_doe_sweep.sh e, para cada
# (kernel, runtime/scheduler), traca GFLOPS medio vs b com barras (min/max),
# DESTACA o b de pico (circulo) e grava a tabela de picos.
#
# A precisao esta embutida no nome do kernel (spotrf/dpotrf/sgetrf_nopiv/
# dgetrf_nopiv), entao facetar por kernel ja separa FP32 de FP64 -- com escala Y
# LIVRE porque FP32 (~15 TFLOPS) e FP64 (~1 TFLOPS) diferem ~10x e no mesmo eixo
# a montanha FP64 vira uma linha reta. Lembrar que N difere por precisao
# (FP32 N=80000, FP64 N=40000), entao NAO comparar b entre precisoes pelo N/b.
#
#   Rscript scripts/analysis/plot_block_size.r [base_dir]
#   (default: data/gpu_tile_poti_* mais recente)
#
# Saida: plots/gflops_vs_b_r.{png,pdf} + plots/block_size_peak.csv

library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) {
  args[1]
} else {
  cands <- sort(Sys.glob("data/gpu_tile_poti_*"), decreasing = TRUE)
  if (length(cands) == 0) {
    stop("nenhum data/gpu_tile_poti_* encontrado; passe base_dir como argumento")
  }
  cands[1]
}
results_file <- file.path(base_dir, "results.csv")
message("lendo ", results_file)

results <- read_csv(results_file, show_col_types = FALSE) |>
  filter(!is.na(gflops))

agg <- results |>
  group_by(runtime, scheduler, algorithm, n, b) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_min  = min(gflops),
    gflops_max  = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = paste(runtime, scheduler, sep = "/"))

# Pico por (kernel, config): o b de maior GFLOPS medio = block size ideal.
peaks <- agg |>
  group_by(algorithm, config, runtime, scheduler, n) |>
  slice_max(gflops_mean, n = 1, with_ties = FALSE) |>
  ungroup()

dir.create("plots", showWarnings = FALSE)
write_csv(
  peaks |>
    select(algorithm, n, runtime, scheduler, b_ideal = b, gflops_peak = gflops_mean) |>
    arrange(algorithm, runtime, scheduler),
  "plots/block_size_peak.csv"
)

node_label <- sub(".*gpu_tile_([a-z]+)_.*", "\\1", base_dir)
gpu_label <- switch(node_label,
  poti = "poti, RTX 4070",
  tupi = "tupi, RTX 4090",
  node_label
)

p <- ggplot(agg, aes(x = b, y = gflops_mean, color = config)) +
  geom_line() +
  geom_point(size = 1) +
  geom_errorbar(
    aes(ymin = gflops_min, ymax = gflops_max),
    width = 0.02 * max(agg$b)
  ) +
  geom_point(data = peaks, size = 3.2, shape = 1, stroke = 1.1) +
  facet_wrap(~algorithm, scales = "free_y", nrow = 2) +
  labs(
    x = "b (tamanho do bloco)", y = "GFLOPS (media; barras = min/max)",
    color = NULL,
    title = sprintf("E1 -- Montanha de block size na GPU (%s)", gpu_label),
    subtitle = "circulo = pico (b ideal). FP32 N=80000, FP64 N=40000"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave("plots/gflops_vs_b_r.png", p, width = 11, height = 7, dpi = 140)
ggsave("plots/gflops_vs_b_r.pdf", p, width = 11, height = 7)
message("escrito plots/gflops_vs_b_r.{png,pdf} e plots/block_size_peak.csv")
