#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "plot_style.r"))

args <- commandArgs(trailingOnly = TRUE)
dirs <- if (length(args) >= 2) {
  args[1:2]
} else if (length(args) == 1) {
  stop("passe 0 ou 2 base_dirs (poti e tupi); recebeu 1")
} else {
  poti_cands <- sort(Sys.glob("data/block_size_poti_*"), decreasing = TRUE)
  tupi_cands <- sort(Sys.glob("data/block_size_tupi_*"), decreasing = TRUE)
  if (length(poti_cands) == 0) stop("nenhum data/block_size_poti_* encontrado")
  if (length(tupi_cands) == 0) stop("nenhum data/block_size_tupi_* encontrado")
  c(poti_cands[1], tupi_cands[1])
}

read_one <- function(base_dir) {
  f <- file.path(base_dir, "results.csv")
  message("lendo ", f)
  node <- basename(base_dir)
  m    <- regexpr("poti|tupi", node)
  if (m > 0) node <- regmatches(node, m)
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

# Teto de b por no: acima disso o overhead nao compensa nessas maquinas.
b_max <- c(poti = 800, tupi = 2000)
results <- results |>
  filter(b <= b_max[as.character(node)])

agg <- results |>
  group_by(node, gpu, runtime, scheduler, algorithm, n, b) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_min  = min(gflops),
    gflops_max  = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = cfg_factor(cfg_label(runtime, scheduler)),
         algo   = algo_label(algorithm))

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

# Largura das barras de erro relativa a grade de b de cada no.
plot_dat <- agg |>
  group_by(node) |>
  mutate(ebw = 0.02 * max(b)) |>
  ungroup()

p <- ggplot(plot_dat, aes(x = b, y = gflops_mean,
                          colour = config, shape = config)) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = gflops_min, ymax = gflops_max, width = ebw)
  ) +
  facet_grid(algo ~ gpu, scales = "free") +
  scale_config_colour() +
  scale_config_shape() +
  expand_y_zero() +
  labs(
    x = "Tamanho do Bloco",
    y = "GFLOPS (média; barras = mín/máx)"
  ) +
  theme_sscad()

ggsave(file.path(out_dir, "gflops_vs_b_compare.png"), p,
       width = FIG_WIDTH_IN, height = 5.9, dpi = 140)
ggsave(file.path(out_dir, "gflops_vs_b_compare.pdf"), p,
       width = FIG_WIDTH_IN, height = 5.9)
