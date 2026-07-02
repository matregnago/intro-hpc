#!/usr/bin/env Rscript
#
# Idea #4 (graficos.md): per-task WAIT TIME -- the gap between a task becoming
# eligible (all its DAG predecessors finished) and actually starting to run:
#
#     wait(T) = StartTime(T) - max(EndTime of T's predecessors)
#
# This is the same ready_time used by parsec_ready_to_parquet.r, kept per-task
# instead of aggregated into the Ready counter. It reads tasks.parquet, whose
# schema is identical for StarPU (native, from tasks.rec) and PaRSEC
# (reconstructed by parsec_tasks_to_parquet.r from .dot + .prof), so ONE code
# path compares both runtimes.
#
# CAVEAT (honest, same as the Ready reconstruction): this is data-availability
# readiness (DAG eligibility), NOT the runtime scheduler's exact queue latency --
# it ignores scheduling overhead. It still measures how long ready work waits
# before a worker picks it up, which is a real scheduler signature.
#
# Outputs:
#   plots/wait_time_dist.{png,pdf}  -- wait distribution per kernel x cfg (violin).
#   plots/wait_time_summary.csv     -- mean/median/p95 wait per cfg x kernel.
#
# Usage:  plot_wait_time.r [base_dir]
#   Processes every run under base_dir that has tasks.parquet. Note: in the
#   `parcial` job only StarPU runs carry tasks.parquet (PaRSEC there lacks the
#   grapher .dot), so it compares StarPU lws vs ws; point it at a job whose
#   PaRSEC runs were reconstructed (e.g. data/manual_traces_*/runs) to get PaRSEC
#   lfq vs gd, or at a job that has both for a cross-runtime comparison.

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

runs <- list_runs(base_dir) %>% filter(.data$tasks)
if (nrow(runs) == 0) {
  stop("no runs with tasks.parquet under ", base_dir,
       " -- for PaRSEC, run parsec_tasks_to_parquet.r first")
}
message("runs com tasks.parquet: ", paste(runs$dir, collapse = ", "))

#' Per-task wait time from a tasks.parquet (DependsOn = space-sep predecessor
#' JobIds). Times normalised to ms (read_tasks_norm) so runtimes are comparable.
wait_times <- function(run_dir) {
  tk <- read_tasks_norm(run_dir) %>%
    filter(!is.na(.data$StartTime), !is.na(.data$EndTime))
  endmap <- setNames(tk$EndTime, tk$JobId)
  gmin   <- min(tk$StartTime)
  ready <- vapply(seq_len(nrow(tk)), function(i) {
    dep <- tk$DependsOn[i]
    if (is.na(dep) || !nzchar(dep)) return(gmin)
    pe <- endmap[strsplit(dep, " ", fixed = TRUE)[[1]]]
    pe <- pe[!is.na(pe)]
    if (!length(pe)) gmin else max(pe)
  }, numeric(1))
  ready <- pmin(ready, tk$StartTime)            # causality clamp
  tk %>%
    mutate(kernel = norm_kernel(.data$Name),
           wait = .data$StartTime - ready,
           start_rel = .data$StartTime - gmin)
}

read_one <- function(row) {
  wait_times(row$path) %>%
    filter(.data$kernel %in% COMPUTE_KERNELS, !.data$kernel %in% INIT_KERNELS) %>%
    transmute(.data$kernel, .data$wait, .data$start_rel, algo = row$algo,
              cfg = paste0(row$runtime, ":", row$scheduler))
}

dat <- bind_rows(lapply(seq_len(nrow(runs)), function(i) read_one(runs[i, ])))
if (nrow(dat) == 0) stop("no compute tasks with wait times")
dat <- dat %>% mutate(panel = paste0(.data$algo, ": ", .data$kernel))

summary_tbl <- dat %>%
  group_by(.data$algo, .data$cfg, .data$kernel) %>%
  summarise(n = n(), mean_wait = mean(.data$wait),
            median_wait = median(.data$wait),
            p95_wait = quantile(.data$wait, 0.95), .groups = "drop") %>%
  arrange(.data$algo, .data$kernel, .data$cfg)
print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(plots_dir(), "wait_time_summary.csv"), row.names = FALSE)

p <- ggplot(dat, aes(.data$cfg, .data$wait + 0.001, fill = .data$cfg)) +
  geom_violin(scale = "width", alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8) +
  scale_y_log10() +
  facet_wrap(~panel, scales = "free_y") +
  labs(title = "Wait time por tarefa (StartTime - prontidao pelo DAG)",
       subtitle = paste0(base_dir, "  |  eixo y = wait+0.001 ms (log)"),
       x = NULL, y = "wait + 0.001 (ms, log10)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
save_plot(p, "wait_time_dist", width = 12, height = 8)
