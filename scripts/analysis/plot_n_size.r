library(ggplot2)
library(dplyr)
library(readr)

# GFLOPS vs N por runtime/scheduler, com faceta por algoritmo (e precisao em runs GPU).
#   Rscript scripts/analysis/plot_n_size.r [results_dir]
# Default: o data/*n_size* mais recente.

args <- commandArgs(trailingOnly = TRUE)
JOB_DIR <- if (length(args) >= 1) {
  args[1]
} else {
  cands <- sort(Sys.glob(c("data/*n_size*", "data/*n*size*")), decreasing = TRUE)
  stopifnot("no data/*n_size* dir found" = length(cands) > 0)
  cands[1]
}
OUT_DIR <- file.path(JOB_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(file.path(JOB_DIR, "results.csv"), show_col_types = FALSE) |>
  filter(!is.na(gflops), gflops > 0)

has_precision <- "precision" %in% names(df)
has_gpus     <- "gpus" %in% names(df) && any(df$gpus > 0, na.rm = TRUE)

df <- df |>
  mutate(
    op        = toupper(sub("^[sdcz]", "", algorithm)),
    config    = paste(runtime, scheduler, sep = "/"),
    precision = if (has_precision) factor(precision, levels = c("FP32", "FP64")) else "FP64"
  )

agg <- df |>
  group_by(runtime, scheduler, algorithm, op, precision, n) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_sd   = sd(gflops),
    n_rep       = dplyr::n(),
    .groups = "drop"
  ) |>
  mutate(
    config    = paste(runtime, scheduler, sep = "/"),
    # configs com 1 rep nao tem DP (NA) -> 0 p/ o ponto ainda plotar sem barra
    gflops_sd = ifelse(is.na(gflops_sd), 0, gflops_sd)
  )

where <- if (has_gpus) {
  sprintf("GPU (%d) | b=%s | t=%s", max(df$gpus, na.rm = TRUE),
          paste(sort(unique(df$b)), collapse = ","),
          paste(sort(unique(df$threads)), collapse = "/"))
} else {
  sprintf("CPU-only | b=%s | t=%s", paste(sort(unique(df$b)), collapse = ","),
          paste(sort(unique(df$threads)), collapse = "/"))
}
ttl <- sprintf("Chameleon %s - GFLOPS vs N", if (has_gpus) "GPU" else "CPU-only")

p <- ggplot(agg, aes(x = n, y = gflops_mean, colour = config)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = gflops_mean - gflops_sd, ymax = gflops_mean + gflops_sd),
                width = 0.03) +
  scale_x_log10() +
  labs(
    x = "N (matriz N x N)", y = "GFLOPS", colour = NULL,
    title = ttl, subtitle = where,
    caption = "Barras de erro: ±1 desvio-padrão entre repetições"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

if (has_precision && nlevels(droplevels(df$precision)) > 1) {
  p <- p + facet_grid(op ~ precision, scales = "free_y")
} else {
  p <- p + facet_wrap(~op, scales = "free_y")
}

out_base <- file.path(OUT_DIR, "gflops_vs_n_r")
ggsave(paste0(out_base, ".png"), p, width = 10, height = 5, dpi = 140)
ggsave(paste0(out_base, ".pdf"), p, width = 10, height = 5)

cat("\nJob:", JOB_DIR, "\n")
cat("Figuras em:", out_base, ".{png,pdf}\n")
