#!/usr/bin/env Rscript
#
# Comparativo poti vs tupi da montanha de block size: GFLOPS medio vs b (barras
# min/max), facet_wrap(algorithm ~ node) com escalas livres, color = config.
# So FP64. As grades de b dos dois jobs nao coincidem (poti: varredura densa;
# tupi: {1000,2000,4000}); no plot a poti corta a cauda b=1000/2000/4000, mas o
# block_size_peak_compare.csv vem da grade completa.
#
#   Rscript scripts/analysis/plot_block_size_compare.r poti_dir tupi_dir
#   (default: data/gpu_tile_poti_* e data/gpu_tile_tupi_* mais recentes)
#
# Saida: plots/final/gflops_vs_b_compare.{png,pdf}
#      + plots/final/block_size_peak_compare.csv

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
dirs <- if (length(args) >= 2) {
  args[1:2]
} else if (length(args) == 1) {
  stop("passe 0 ou 2 base_dirs (poti e tupi); recebeu 1")
} else {
  poti_cands <- sort(Sys.glob("data/gpu_tile_poti_*"), decreasing = TRUE)
  tupi_cands <- sort(Sys.glob("data/gpu_tile_tupi_*"), decreasing = TRUE)
  if (length(poti_cands) == 0) stop("nenhum data/gpu_tile_poti_* encontrado")
  if (length(tupi_cands) == 0) stop("nenhum data/gpu_tile_tupi_* encontrado")
  c(poti_cands[1], tupi_cands[1])
}

read_one <- function(base_dir) {
  f <- file.path(base_dir, "results.csv")
  message("lendo ", f)
  node <- sub(".*gpu_tile_([a-z]+)_.*", "\\1", base_dir)
  gpu <- switch(node,
    poti = "poti (RTX 4070)",
    tupi = "tupi (RTX 4090)",
    node)
  read_csv(f, show_col_types = FALSE) |>
    filter(!is.na(gflops)) |>
    mutate(node = factor(node, levels = c("poti", "tupi")),
           gpu  = gpu)
}

results <- bind_rows(lapply(dirs, read_one)) |>
  filter(precision == "FP64")

agg <- results |>
  group_by(node, gpu, runtime, scheduler, algorithm, n, b) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_min  = min(gflops),
    gflops_max  = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = paste(runtime, scheduler, sep = "/"))

# Pico por (node, kernel, config): o b de maior GFLOPS medio = b ideal.
peaks <- agg |>
  group_by(node, gpu, algorithm, config, runtime, scheduler, n) |>
  slice_max(gflops_mean, n = 1, with_ties = FALSE) |>
  ungroup()

out_dir <- "plots/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(
  peaks |>
    select(node, gpu, algorithm, n, runtime, scheduler,
           b_ideal = b, gflops_peak = gflops_mean) |>
    arrange(node, algorithm, runtime, scheduler),
  file.path(out_dir, "block_size_peak_compare.csv")
)

# So no plot: poti sem a cauda plana, largura das barras relativa a grade de b
# de cada no (as grades diferem ~4x).
plot_dat <- agg |>
  filter(!(node == "poti" & b %in% c(1000, 2000, 4000))) |>
  group_by(node) |>
  mutate(ebw = 0.02 * max(b)) |>
  ungroup()

# scales="free": as grades de b e o teto de GFLOPS (4070 capada em FP64 vs
# 4090) diferem demais para compartilhar eixos.
p <- ggplot(plot_dat, aes(x = b, y = gflops_mean, color = config)) +
  geom_line() +
  geom_point(size = 1) +
  geom_errorbar(
    aes(ymin = gflops_min, ymax = gflops_max, width = ebw)
  ) +
  facet_wrap(vars(algorithm, gpu), scales = "free", ncol = 2) +
  labs(
    x = "b (tamanho do bloco)",
    y = "GFLOPS (media; barras = min/max)",
    color = NULL,
    subtitle = "FP64, N=40000"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "gflops_vs_b_compare.png"), p,
       width = 13, height = 7, dpi = 140)
ggsave(file.path(out_dir, "gflops_vs_b_compare.pdf"), p, width = 13, height = 7)
