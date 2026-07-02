#!/usr/bin/env Rscript
#
# Idea #1 (graficos.md) + extra: per-task-type timing, compared across runtimes
# and schedulers. Reads each run's application.parquet, keeps the compute kernels
# (normalised so StarPU's dgemm/dpotrf line up with PaRSEC's gemm/potrf), and
# renders two views:
#
#   plots/task_times_mean.{png,pdf}    -- mean kernel duration, bars grouped by
#                                         runtime:scheduler, facet per kernel.
#   plots/task_times_violin.{png,pdf}  -- full duration distribution per kernel
#                                         (exposes jitter/variance, e.g. a kernel
#                                         that serialises in one runtime).
# Also prints a summary table (mean/median/sd/n per runtime:scheduler x kernel).
#
# Source per run is runtime-agnostic via read_exec(): application.parquet where
# available (StarPU; the CPU PaRSEC job), else the GPU-complete tasks.parquet
# (PaRSEC GPU). Faceted by algorithm, so a job mixing dpotrf/dgeqrf (or the GPU
# spotrf/sgeqrf job) plots cleanly.
#
# On GPU jobs the same kernel runs in two very different modes, so everything is
# additionally split by resource class (CPU vs GPU, from the ResourceId/worker of
# each span). Aggregating the two classes hides the finding that matters on the
# poti FP64 job: StarPU's CPU gemms are ~25% slower than PaRSEC's (its 4 dedicated
# CUDA-worker threads oversubscribe the cores), which is where dmda loses the
# makespan it gains on the GPU. NB: PaRSEC's GPU spans are host-observed
# (submit->completion event) and overlap on the 4 streams, so GPU-class durations
# are inflated for PaRSEC; the CPU-vs-CPU comparison is the like-for-like one.
#
# Usage:
#   plot_task_times.r [base_dir] [algo]
# base_dir defaults to the GPU job (spotrf/sgeqrf n19200 b960, 8 runs).
# algo (optional) keeps only that algorithm's runs (e.g. dpotrf); with a single
# algorithm the facet labels drop the "algo: " prefix.

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(forcats)
})

# locate this script's dir to source the shared helpers
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"
algo_sel <- if (length(args) >= 2) args[[2]] else NA

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (!is.na(algo_sel)) runs <- runs %>% filter(.data$algo == algo_sel)
if (nrow(runs) == 0) stop("no runs with application/tasks parquet under ", base_dir,
                          if (!is.na(algo_sel)) paste0(" (algo=", algo_sel, ")") else "")
message("runs: ", paste(runs$dir, collapse = ", "))

# Per run: unified compute-kernel intervals, tagged with experiment metadata.
read_kernels <- function(row) {
  ex <- read_exec(row$path)
  if (is.null(ex)) return(NULL)
  ex %>%
    transmute(.data$kernel, .data$Duration, .data$src,
              class = ifelse(grepl("CUDA", .data$worker), "GPU", "CPU"),
              runtime = row$runtime, scheduler = row$scheduler,
              algo = row$algo, cfg = paste0(row$runtime, ":", row$scheduler))
}

dat <- bind_rows(lapply(seq_len(nrow(runs)), function(i) read_kernels(runs[i, ])))
if (nrow(dat) == 0) stop("no compute-kernel states found")
message("fontes: ", paste(unique(paste0(dat$runtime, "=", dat$src)), collapse = ", "))

# drop init/util kernels from the headline view -- the factorization kernels are
# what Schnorr's "mean time per task type" is about; keep them in the summary.
core <- dat %>% filter(!.data$kernel %in% INIT_KERNELS)

summary_tbl <- dat %>%
  group_by(.data$algo, .data$cfg, .data$class, .data$kernel) %>%
  summarise(n = n(), mean_ms = mean(.data$Duration),
            median_ms = median(.data$Duration), sd_ms = sd(.data$Duration),
            .groups = "drop") %>%
  arrange(.data$algo, .data$kernel, .data$cfg, .data$class)
print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(plots_dir(), "task_times_summary.csv"), row.names = FALSE)

# facet label combines algorithm + kernel (kernels differ per algorithm); with a
# single algorithm the prefix is redundant, so it's just the kernel.
one_algo <- n_distinct(core$algo) == 1
core <- core %>%
  mutate(panel = if (one_algo) .data$kernel
         else paste0(.data$algo, ": ", .data$kernel))

caption_cls <- paste("classe = onde o span executou (ResourceId).",
                     "GPU do PaRSEC: spans host-observados com overlap entre streams",
                     "(duracao inflada) -- comparar GPU só entre configs StarPU;",
                     "CPU vs CPU e like-for-like.")

# ---- (1) mean duration bars, facet per algo+kernel, split CPU/GPU ----
mean_tbl <- core %>%
  group_by(.data$panel, .data$cfg, .data$class) %>%
  summarise(mean_us = mean(.data$Duration),
            se = sd(.data$Duration) / sqrt(n()), .groups = "drop")

p_mean <- ggplot(mean_tbl, aes(.data$cfg, .data$mean_us, fill = .data$class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(aes(ymin = .data$mean_us - .data$se,
                    ymax = .data$mean_us + .data$se),
                position = position_dodge(width = 0.8), width = 0.3) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c(CPU = "#377eb8", GPU = "#e41a1c")) +
  labs(title = "Tempo medio por tipo de tarefa (kernel) x runtime:scheduler x classe",
       x = NULL, y = "duracao media (ms)", fill = "classe") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")
save_plot(p_mean, "task_times_mean", width = 12, height = 8)

# ---- (2) duration distribution (violin + boxplot) per algo+kernel ----
p_violin <- ggplot(core, aes(.data$cfg, .data$Duration, fill = .data$class)) +
  geom_violin(scale = "width", alpha = 0.5, colour = NA,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c(CPU = "#377eb8", GPU = "#e41a1c")) +
  labs(title = "Distribuicao da duracao por tarefa (kernel) x runtime:scheduler x classe",
       caption = caption_cls,
       x = NULL, y = "duracao (ms)", fill = "classe") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")
save_plot(p_violin, "task_times_violin", width = 12, height = 8)
