library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# E2 -- Impacto do escalonador na GPU: StarPU vs PaRSEC.
# Analogo ao distrib-scheduler.pdf do Schnorr (WAMCA 2025): violino de makespan
# por escalonador, com runtime na cor e algoritmo x precisao nas facetas.
# Tambem emite a versao em GFLOPS e uma tabela-resumo (media +/- dp + cv).
#
#   Rscript scripts/analysis/plot_gpu_sched.r [results_dir]
# Default: o data/gpu_sched_* mais recente.

args <- commandArgs(trailingOnly = TRUE)
JOB_DIR <- if (length(args) >= 1) {
  args[1]
} else {
  cands <- sort(Sys.glob("data/gpu_sched_*"), decreasing = TRUE)
  stopifnot("no data/gpu_sched_* dir found" = length(cands) > 0)
  cands[1]
}
OUT_DIR <- file.path(JOB_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(file.path(JOB_DIR, "results.csv"), show_col_types = FALSE) |>
  filter(!is.na(gflops), gflops > 0) |>
  mutate(
    op        = toupper(sub("^[sdcz]", "", algorithm)),
    precision = factor(precision, levels = c("FP32", "FP64")),
    runtime   = factor(runtime),
    scheduler = factor(scheduler)
  )

sub_lab <- sprintf("GPU on | N=%s | b=%s | t=%s | reps=%d",
                   paste(sort(unique(df$n)), collapse = ","),
                   paste(sort(unique(df$b)), collapse = ","),
                   paste(sort(unique(df$threads)), collapse = "/"),
                   max(df$rep))

violin <- function(yvar, ylab, file, title) {
  p <- ggplot(df, aes(x = scheduler, y = .data[[yvar]], fill = runtime)) +
    geom_violin(scale = "width", alpha = 0.55, colour = NA,
                position = position_dodge(width = 0.8)) +
    geom_jitter(aes(colour = runtime), width = 0.12, size = 0.6, alpha = 0.5,
                show.legend = FALSE) +
    facet_grid(op ~ precision, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Escalonador", y = ylab, fill = NULL, title = title,
         subtitle = sub_lab) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUT_DIR, paste0(file, ".png")), p, width = 9, height = 6, dpi = 140)
  ggsave(file.path(OUT_DIR, paste0(file, ".pdf")), p, width = 9, height = 6)
}

violin("time", "Makespan [s]", "makespan_by_scheduler",
       "Impacto do escalonador na GPU: StarPU vs PaRSEC")
violin("gflops", "GFLOPS", "gflops_by_scheduler",
       "Desempenho por escalonador na GPU: StarPU vs PaRSEC")

# --- Tabela-resumo: media +/- dp + coeficiente de variacao ------------------
summary_tbl <- df |>
  group_by(precision, op, runtime, scheduler) |>
  summarise(
    n_runs      = dplyr::n(),
    time_mean   = round(mean(time), 3),
    time_sd     = round(sd(time), 3),
    gflops_mean = round(mean(gflops), 1),
    gflops_sd   = round(sd(gflops), 1),
    cv_pct      = round(100 * sd(time) / mean(time), 1),
    .groups = "drop"
  ) |>
  arrange(precision, op, desc(gflops_mean))

write_csv(summary_tbl, file.path(OUT_DIR, "summary_gpu_sched.csv"))

# --- Melhor escalonador por celula (precision x op x runtime) ---------------
best_tbl <- summary_tbl |>
  group_by(precision, op, runtime) |>
  slice_max(gflops_mean, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(precision, op, runtime, best_sched = scheduler, gflops_mean)

write_csv(best_tbl, file.path(OUT_DIR, "best_scheduler.csv"))

cat("\n=== Resumo (media sobre reps) ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\n=== Melhor escalonador por (precisao, algoritmo, runtime) ===\n")
print(as.data.frame(best_tbl), row.names = FALSE)
cat("\nFiguras + tabelas em:", OUT_DIR, "\n")
