#!/usr/bin/env Rscript
#
# Per-task-type timing compared across runtimes/schedulers, from each run's
# application.parquet (fallback: tasks.parquet, the GPU-complete PaRSEC source
# -- see read_exec). Splits every kernel by resource class (CPU vs GPU): NB
# PaRSEC's GPU spans are host-observed and overlap across streams, so
# GPU-class durations are inflated for PaRSEC; CPU vs CPU is like-for-like.
#
# Usage:  plot_task_times.r [base_dir] [algo]
# Output: plots/task_times_mean.{png,pdf}, plots/task_times_violin.{png,pdf}
#         + task_times_summary.csv. Honours PLOTS_DIR.

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(forcats)
})

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

# init/util kernels stay in the summary but out of the plots.
core <- dat %>% filter(!.data$kernel %in% INIT_KERNELS)

summary_tbl <- dat %>%
  group_by(.data$algo, .data$cfg, .data$class, .data$kernel) %>%
  summarise(n = n(), mean_ms = mean(.data$Duration),
            median_ms = median(.data$Duration), sd_ms = sd(.data$Duration),
            .groups = "drop") %>%
  arrange(.data$algo, .data$kernel, .data$cfg, .data$class)
print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(plots_dir(), "task_times_summary.csv"), row.names = FALSE)

# facet label = algo + kernel; with a single algorithm just the kernel.
one_algo <- n_distinct(core$algo) == 1
core <- core %>%
  mutate(panel = if (one_algo) .data$kernel
         else paste0(.data$algo, ": ", .data$kernel))

caption_cls <- paste("classe = onde o span executou (ResourceId).",
                     "GPU do PaRSEC: spans host-observados com overlap entre streams",
                     "(duracao inflada) -- comparar GPU só entre configs StarPU;",
                     "CPU vs CPU e like-for-like.")

# (1) mean duration bars, facet per algo+kernel, split CPU/GPU
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

# (2) duration distribution (violin + boxplot) per algo+kernel
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
