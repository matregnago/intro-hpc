library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# CPU-vs-GPU validation job (slurm/cuda_vs_cpu.slurm). Single-precision sweep
# (sgemm/spotrf/sgeqrf) on the poti RTX 4070, one scheduler per runtime
# (starpu:lws, parsec:lfq), modes cpu (16 threads, 0 gpu) and gpu (8 threads,
# 1 gpu). Edit JOB_DIR to point at a new job id.
JOB_DIR <- "data/chameleon-cuda-vscpu_791531"
OUT_DIR <- file.path(JOB_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

results <- read_csv(file.path(JOB_DIR, "results.csv"), show_col_types = FALSE)

# Rótulos derivados dos dados, p/ o script servir a qualquer job desta família.
b_lab   <- paste(sort(unique(results$b)), collapse = ",")
cpu_thr <- unique(results$threads[results$mode == "cpu"])
gpu_thr <- unique(results$threads[results$mode == "gpu"])
sub_lab <- sprintf("cpu: %s threads / 0 gpu    gpu: %s threads / 1 gpu",
                   paste(cpu_thr, collapse = "/"), paste(gpu_thr, collapse = "/"))

# Mean / spread per cell. config = runtime + mode is the comparison line.
agg <- results |>
  group_by(runtime, mode, scheduler, algorithm, n) |>
  summarise(
    time_mean   = mean(time),
    gflops_mean = mean(gflops),
    gflops_min  = min(gflops),
    gflops_max  = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = paste(runtime, mode, sep = "/"))

# --- Plot 1: GFLOPS vs N, cpu vs gpu, faceted by algorithm -----------------
p_gflops <- ggplot(agg, aes(x = n, y = gflops_mean, color = config)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = gflops_min, ymax = gflops_max),
                width = 0.02 * max(agg$n)) +
  facet_wrap(~algorithm, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(agg$n))) +
  labs(
    x = "N (matriz N x N)", y = "GFLOPS", color = NULL,
    title = sprintf("Chameleon CPU vs GPU (poti, RTX 4070) - FP32, b=%s", b_lab),
    subtitle = sub_lab
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "gflops_vs_n_cuda_vscpu.png"), p_gflops,
       width = 11, height = 4.2, dpi = 140)
ggsave(file.path(OUT_DIR, "gflops_vs_n_cuda_vscpu.pdf"), p_gflops,
       width = 11, height = 4.2)

# --- Speedup: gpu vs cpu within each runtime -------------------------------
speedup <- agg |>
  select(runtime, mode, algorithm, n, gflops_mean, time_mean) |>
  pivot_wider(names_from = mode,
              values_from = c(gflops_mean, time_mean)) |>
  mutate(speedup = gflops_mean_gpu / gflops_mean_cpu)

p_speedup <- ggplot(speedup,
                    aes(x = factor(n), y = speedup, fill = runtime)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.1fx", speedup)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  facet_wrap(~algorithm) +
  labs(
    x = "N (matriz N x N)", y = "Speedup GPU / CPU", fill = NULL,
    title = "Aceleracao GPU vs CPU (poti, RTX 4070) - FP32",
    subtitle = "linha tracejada = sem ganho (1x)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "speedup_gpu_vs_cpu.png"), p_speedup,
       width = 11, height = 4.2, dpi = 140)
ggsave(file.path(OUT_DIR, "speedup_gpu_vs_cpu.pdf"), p_speedup,
       width = 11, height = 4.2)

# --- Tidy summary table ----------------------------------------------------
summary_tbl <- speedup |>
  transmute(
    runtime, algorithm, n,
    gflops_cpu = round(gflops_mean_cpu, 1),
    gflops_gpu = round(gflops_mean_gpu, 1),
    time_cpu_s = round(time_mean_cpu, 3),
    time_gpu_s = round(time_mean_gpu, 3),
    speedup    = round(speedup, 2)
  ) |>
  arrange(algorithm, runtime, n)

write_csv(summary_tbl, file.path(OUT_DIR, "summary_cuda_vscpu.csv"))

cat("\n=== CPU vs GPU summary (mean over reps) ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\nWrote:\n")
cat(" ", file.path(OUT_DIR, "gflops_vs_n_cuda_vscpu.{png,pdf}"), "\n")
cat(" ", file.path(OUT_DIR, "speedup_gpu_vs_cpu.{png,pdf}"), "\n")
cat(" ", file.path(OUT_DIR, "summary_cuda_vscpu.csv"), "\n")
