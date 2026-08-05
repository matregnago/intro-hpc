#!/usr/bin/env Rscript
#
# Comparativo poti vs tupi de GFLOPS vs N (media +- 1 desvio-padrao entre
# repeticoes), facet_grid(op ~ node), color = config. Cada no roda no seu b
# ideal; o subtitulo deriva isso dos dados.
#
#   Rscript scripts/analysis/plot_n_size_compare.r poti_dir tupi_dir
#   (default: data/n_size_poti_* e data/n_size_tupi_* mais recentes)
#
# Saida: plots/final/gflops_vs_n_compare.{png,pdf}

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
  poti_cands <- sort(Sys.glob("data/n_size_poti_*"), decreasing = TRUE)
  tupi_cands <- sort(Sys.glob("data/n_size_tupi_*"), decreasing = TRUE)
  if (length(poti_cands) == 0) stop("nenhum data/n_size_poti_* encontrado")
  if (length(tupi_cands) == 0) stop("nenhum data/n_size_tupi_* encontrado")
  c(poti_cands[1], tupi_cands[1])
}

read_one <- function(base_dir) {
  f <- file.path(base_dir, "results.csv")
  message("lendo ", f)
  # node e extraido do nome do dir independente do prefixo; se nao casa, usa o
  # proprio nome.
  node <- basename(base_dir)
  m    <- regexpr("poti|tupi", node)
  if (m > 0) node <- regmatches(node, m)
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
  mutate(op = algo_label(algorithm))

agg <- results |>
  group_by(node, gpu, runtime, scheduler, op, n) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_sd   = sd(gflops),
    .groups = "drop"
  ) |>
  mutate(
    config    = cfg_factor(cfg_label(runtime, scheduler)),
    # 1 rep -> sd NA -> 0, p/ o ponto plotar sem barra
    gflops_sd = ifelse(is.na(gflops_sd), 0, gflops_sd)
  ) |>
  # Largura da barra de erro em unidades de N: com o eixo X linear, uma largura
  # fixa some. 1.5% da amplitude de N do proprio no.
  group_by(node) |>
  mutate(ebw = 0.015 * diff(range(n))) |>
  ungroup()

# b ideal de cada no (constante dentro do job). Nao vira subtitulo: o revisor
# pediu que o conteudo do titulo va para o caption do artigo -- imprimimos aqui
# para que o texto do caption possa ser escrito a partir dos dados.
b_tbl <- results |>
  distinct(node, b) |>
  arrange(node)
message("para o caption: FP64, ",
        paste(sprintf("%s b=%d", b_tbl$node, b_tbl$b), collapse = ", "))

# scales="free_y" libera o eixo por linha; as duas colunas (poti | tupi)
# compartilham o eixo, entao a comparacao entre nos fica direta.
# Eixo X LINEAR (era log10): em log a distancia entre os niveis de N engana a
# leitura da curva. Eixo Y inclui o zero.
p <- ggplot(agg, aes(x = n, y = gflops_mean,
                     colour = config, shape = config)) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = gflops_mean - gflops_sd,
                    ymax = gflops_mean + gflops_sd,
                    width = ebw)) +
  facet_grid(op ~ gpu, scales = "free_y") +
  scale_config_colour() +
  scale_config_shape() +
  expand_y_zero() +
  labs(
    x = "N (matriz N x N)",
    y = "GFLOPS (média)"
  ) +
  theme_sscad()

out_dir <- "plots/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "gflops_vs_n_compare.png"), p,
       width = 13, height = 7, dpi = 140)
ggsave(file.path(out_dir, "gflops_vs_n_compare.pdf"), p, width = 13, height = 7)
message("escrito ", out_dir, "/gflops_vs_n_compare.{png,pdf}")
