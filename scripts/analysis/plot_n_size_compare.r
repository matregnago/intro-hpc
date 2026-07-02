#!/usr/bin/env Rscript

# Comparativo poti vs tupi de GFLOPS vs N. Mesma metrica do plot_n_size.r
# (media +- 1 desvio-padrao entre repeticoes), mas sobrepoe os dois nos no mesmo
# painel -- facet_grid(op ~ node), color = config. Cada no roda no seu b ideal
# (E1): poti b=500, tupi b=1000 -- o subtitulo deriva isso dos dados.
#
#   Rscript scripts/analysis/plot_n_size_compare.r poti_dir tupi_dir
#   (default: data/gpu_n_size_poti_* e data/gpu_n_size_tupi_* mais recentes)
#
# Saida: plots/final/gflops_vs_n_compare.{png,pdf}

library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
dirs <- if (length(args) >= 2) {
  args[1:2]
} else if (length(args) == 1) {
  stop("passe 0 ou 2 base_dirs (poti e tupi); recebeu 1")
} else {
  poti_cands <- sort(Sys.glob("data/gpu_n_size_poti_*"), decreasing = TRUE)
  tupi_cands <- sort(Sys.glob("data/gpu_n_size_tupi_*"), decreasing = TRUE)
  if (length(poti_cands) == 0) stop("nenhum data/gpu_n_size_poti_* encontrado")
  if (length(tupi_cands) == 0) stop("nenhum data/gpu_n_size_tupi_* encontrado")
  c(poti_cands[1], tupi_cands[1])
}

read_one <- function(base_dir) {
  f <- file.path(base_dir, "results.csv")
  message("lendo ", f)
  node <- sub(".*gpu_n_size_([a-z]+)_.*", "\\1", base_dir)
  gpu <- switch(node,
    poti = "poti (RTX 4070)",
    tupi = "tupi (RTX 4090)",
    node)
  read_csv(f, show_col_types = FALSE) |>
    filter(!is.na(gflops), gflops > 0) |>
    mutate(node = factor(node, levels = c("poti", "tupi")),
           gpu  = gpu)
}

results <- bind_rows(lapply(dirs, read_one)) |>
  filter(precision == "FP64") |>
  mutate(op = toupper(sub("^[sdcz]", "", algorithm)))

agg <- results |>
  group_by(node, gpu, runtime, scheduler, op, n) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_sd   = sd(gflops),
    .groups = "drop"
  ) |>
  mutate(
    config    = paste(runtime, scheduler, sep = "/"),
    # configs com 1 rep nao tem DP (NA) -> 0 p/ o ponto ainda plotar sem barra
    gflops_sd = ifelse(is.na(gflops_sd), 0, gflops_sd)
  )

# b ideal de cada no (constante dentro do job), p/ o subtitulo
b_tbl <- results |>
  distinct(node, b) |>
  arrange(node)
sub <- paste0("FP64  |  ",
              paste(sprintf("%s: b=%d", b_tbl$node, b_tbl$b), collapse = "  |  "))

# facet_grid(op ~ gpu): colunas = nos (poti | tupi), linhas = algoritmos.
# scales="free_y" libera o eixo por LINHA; cada linha compartilha Y entre as
# duas colunas, entao a comparacao poti-vs-tupi fica no mesmo eixo.
p <- ggplot(agg, aes(x = n, y = gflops_mean, colour = config)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = gflops_mean - gflops_sd, ymax = gflops_mean + gflops_sd),
                width = 0.03) +
  scale_x_log10() +
  facet_grid(op ~ gpu, scales = "free_y") +
  labs(
    x = "N (matriz N x N)",
    y = "GFLOPS (media)",
    colour = NULL,
    subtitle = sub
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

out_dir <- "plots/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "gflops_vs_n_compare.png"), p,
       width = 13, height = 7, dpi = 140)
ggsave(file.path(out_dir, "gflops_vs_n_compare.pdf"), p, width = 13, height = 7)
message("escrito ", out_dir, "/gflops_vs_n_compare.{png,pdf}")
