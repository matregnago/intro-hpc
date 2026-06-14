library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# StarPU-vs-PaRSEC GPU comparison (scripts/run_local_gpu_compare.sh).
# Both runtimes, GPU enabled, {dpotrf,dgeqrf} over a sweep of N, double
# precision. Usage:
#   Rscript scripts/analysis/plot_local_gpu_compare.r <results_dir>
# Defaults to the most recent data/local_gpu_compare_* if no arg is given.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  JOB_DIR <- args[1]
} else {
  cands <- sort(Sys.glob("data/local_gpu_compare_*"), decreasing = TRUE)
  stopifnot("no data/local_gpu_compare_* dir found" = length(cands) > 0)
  JOB_DIR <- cands[1]
}
OUT_DIR <- file.path(JOB_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

results <- read_csv(file.path(JOB_DIR, "results.csv"), show_col_types = FALSE) |>
  filter(!is.na(gflops))

# Precision is encoded in the algorithm prefix: s=FP32, d=FP64, c/z=complex.
prec_of <- function(algo) {
  p <- substr(algo, 1, 1)
  c(s = "FP32", d = "FP64", c = "complex FP32", z = "complex FP64")[p] |> unname()
}
prec_lab <- paste(sort(unique(prec_of(results$algorithm))), collapse = "/")

b_lab   <- paste(sort(unique(results$b)), collapse = ",")
thr_lab <- paste(sort(unique(results$threads)), collapse = "/")
gpu_lab <- paste(sort(unique(results$gpus)), collapse = "/")
sub_lab <- sprintf("GPU habilitado | %s | b=%s | t=%s | g=%s | reps=%d (RTX 4060 Ti, i5-14600K)",
                   prec_lab, b_lab, thr_lab, gpu_lab, max(results$rep))

# Mean / spread per cell. The comparison line is the runtime.
agg <- results |>
  group_by(runtime, scheduler, algorithm, n) |>
  summarise(
    time_mean   = mean(time),
    gflops_mean = mean(gflops),
    gflops_min  = min(gflops),
    gflops_max  = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = sprintf("%s (%s)", runtime, scheduler))

# --- Plot 1: GFLOPS vs N, starpu vs parsec, faceted by algorithm -----------
p_gflops <- ggplot(agg, aes(x = n, y = gflops_mean, color = config)) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = gflops_min, ymax = gflops_max),
                width = 0.02 * max(agg$n)) +
  facet_wrap(~algorithm, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(agg$n))) +
  labs(
    x = "N (matriz N x N)", y = "GFLOPS", color = NULL,
    title = sprintf("Chameleon GPU: StarPU vs PaRSEC (%s)", prec_lab),
    subtitle = sub_lab
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "gflops_vs_n_starpu_vs_parsec.png"), p_gflops,
       width = 11, height = 4.4, dpi = 140)
ggsave(file.path(OUT_DIR, "gflops_vs_n_starpu_vs_parsec.pdf"), p_gflops,
       width = 11, height = 4.4)

# --- Plot 2: PaRSEC / StarPU GFLOPS ratio ----------------------------------
ratio <- agg |>
  select(runtime, algorithm, n, gflops_mean) |>
  pivot_wider(names_from = runtime, values_from = gflops_mean)

if (all(c("starpu", "parsec") %in% names(ratio))) {
  ratio <- ratio |> mutate(ratio = parsec / starpu)

  p_ratio <- ggplot(ratio, aes(x = factor(n), y = ratio)) +
    geom_col(fill = "#3b6ea5", width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = sprintf("%.2fx", ratio)), vjust = -0.3, size = 3) +
    facet_wrap(~algorithm) +
    labs(
      x = "N (matriz N x N)", y = "GFLOPS PaRSEC / StarPU",
      title = sprintf("Razao de desempenho PaRSEC vs StarPU na GPU (%s)", prec_lab),
      subtitle = "linha tracejada = empate (1x); >1 favorece PaRSEC"
    ) +
    theme_bw()

  ggsave(file.path(OUT_DIR, "ratio_parsec_vs_starpu.png"), p_ratio,
         width = 11, height = 4.2, dpi = 140)
  ggsave(file.path(OUT_DIR, "ratio_parsec_vs_starpu.pdf"), p_ratio,
         width = 11, height = 4.2)
}

# --- Tidy summary table ----------------------------------------------------
summary_tbl <- agg |>
  transmute(runtime, scheduler, algorithm, n,
            gflops = round(gflops_mean, 1),
            time_s = round(time_mean, 3)) |>
  arrange(algorithm, runtime, n)

write_csv(summary_tbl, file.path(OUT_DIR, "summary_local_gpu_compare.csv"))

cat("\n=== StarPU vs PaRSEC GPU summary (mean over reps) ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\nWrote plots + summary to:", OUT_DIR, "\n")
