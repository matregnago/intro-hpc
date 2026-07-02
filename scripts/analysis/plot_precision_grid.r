#!/usr/bin/env Rscript
#
# A contraprova do "PaRSEC lida melhor com GPU" em uma figura: o melhor GFLOPS
# medio de cada runtime:scheduler, por algoritmo, em tres regimes de hardware --
#
#   poti (RTX 4070)  FP64: pico FP64 da GPU (458 GF, crippled 1:64) MENOR que o
#                          agregado das CPUs -> unico regime onde parsec:gd lidera;
#   tupi (RTX 4090)  FP64: pico FP64 ~1.29 TF ja MAIOR que as CPUs -> StarPU lidera;
#   poti (RTX 4070)  FP32: pico 29.3 TF, GPU domina -> StarPU dmda/dmdas 2.3x acima.
#
# Ou seja: o gd so vence onde a jogada certa e quase ignorar a GPU. Quando a GPU
# vale a pena, o perf-model do dmda faz diferenca e o device layer estatico do
# PaRSEC fica para tras.
#
# Para cada (maquina, precisao, algoritmo, config) toma-se a MELHOR celula do
# sweep: media das repeticoes por (n, b), depois o maximo sobre as celulas --
# "o melhor que aquela config consegue naquele regime", imune a celulas
# desfavoraveis que so um dos sweeps visitou. Linha tracejada = pico FP64/FP32
# da GPU (specs; FP64 = FP32/64 nas GeForce Ada).
#
# Uso:  plot_precision_grid.r [job_dir ...]
# Default: os 4 jobs gpu_{n_size,tile}_{poti,tupi} mais recentes em data/.

suppressMessages({
  library(dplyr); library(ggplot2); library(readr); library(stringr); library(tidyr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args <- commandArgs(trailingOnly = TRUE)
jobs <- if (length(args) >= 1) args else {
  pick_latest <- function(glob) {
    cands <- sort(Sys.glob(glob), decreasing = TRUE)
    if (length(cands)) cands[1] else NULL
  }
  unlist(lapply(c("data/gpu_n_size_poti_*", "data/gpu_tile_poti_*",
                  "data/gpu_n_size_tupi_*", "data/gpu_tile_tupi_*"), pick_latest))
}
stopifnot("no result jobs found" = length(jobs) > 0)
message("jobs: ", paste(basename(jobs), collapse = ", "))

GPU_META <- tribble(
  ~machine, ~gpu,       ~fp64_peak, ~fp32_peak,
  "poti",   "RTX 4070",        458,      29300,
  "tupi",   "RTX 4090",       1290,      82600
)

df <- bind_rows(lapply(jobs, function(j) {
  read_csv(file.path(j, "results.csv"), show_col_types = FALSE) %>%
    mutate(machine = str_match(basename(j), "(poti|tupi)")[, 2])
})) %>%
  filter(!is.na(.data$gflops)) %>%
  left_join(GPU_META, by = "machine") %>%
  mutate(
    algo_stem = sub("^[sd]", "", .data$algorithm),
    cfg       = paste0(.data$runtime, ":", .data$scheduler),
    gpu_peak  = ifelse(.data$precision == "FP64", .data$fp64_peak, .data$fp32_peak),
    regime    = sprintf("%s (%s) - %s | pico GPU %s",
                        .data$machine, .data$gpu, .data$precision,
                        ifelse(.data$gpu_peak >= 1000,
                               sprintf("%.1f TF", .data$gpu_peak / 1000),
                               sprintf("%d GF", .data$gpu_peak)))
  )

# ordena os regimes por pico de GPU crescente: a historia e "conforme a GPU
# passa a valer a pena, o StarPU assume a lideranca"
regime_lvls <- df %>% distinct(.data$regime, .data$gpu_peak) %>%
  arrange(.data$gpu_peak) %>% pull(.data$regime)
df <- df %>% mutate(regime = factor(.data$regime, levels = regime_lvls))

best <- df %>%
  group_by(.data$machine, .data$regime, .data$gpu_peak, .data$algo_stem,
           .data$runtime, .data$cfg, .data$n, .data$b) %>%
  summarise(gf = mean(.data$gflops), .groups = "drop") %>%
  group_by(.data$machine, .data$regime, .data$gpu_peak, .data$algo_stem,
           .data$runtime, .data$cfg) %>%
  slice_max(.data$gf, n = 1, with_ties = FALSE) %>%
  ungroup()

print(as.data.frame(best %>% arrange(.data$regime, .data$algo_stem, desc(.data$gf))),
      digits = 5)
write.csv(best, file.path(plots_dir(), "precision_grid_summary.csv"),
          row.names = FALSE)

peaks <- best %>% distinct(.data$regime, .data$algo_stem, .data$gpu_peak)

p <- ggplot(best, aes(.data$cfg, .data$gf, fill = .data$runtime)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_hline(data = peaks, aes(yintercept = .data$gpu_peak),
             linetype = "dashed", colour = "grey30") +
  geom_text(aes(label = sprintf("%.0f", .data$gf)), vjust = -0.35, size = 3.1) +
  geom_text(aes(y = 0, label = sprintf("n=%d b=%d", .data$n, .data$b)),
            vjust = 1.4, size = 2.4, colour = "grey30") +
  facet_wrap(~algo_stem + regime, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(parsec = "#E69F00", starpu = "#56B4E9")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  coord_cartesian(clip = "off") +
  labs(title = "Melhor GFLOPS por config em tres regimes de GPU",
       subtitle = paste("melhor media do sweep (n, b) por config;",
                        "tracejado = pico teorico da GPU na precisao do painel"),
       caption = paste("parsec:gd so lidera no regime em que o pico FP64 da GPU",
                       "(458 GF) fica abaixo do agregado das CPUs -- onde a",
                       "estrategia otima e quase nao usar a GPU"),
       x = NULL, y = "GFLOPS (melhor media)", fill = "runtime") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom")
save_plot(p, "precision_grid", width = 13, height = 8)
