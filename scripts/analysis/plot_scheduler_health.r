#!/usr/bin/env Rscript
#
# Concurrency profile to differentiate the schedulers: number of WORKERS busy with
# an application factorization kernel over time (FACTORIZATION_KERNELS only;
# matrix-gen/init dropped). Shows ramp-up / steady-state / tail and how well each
# runtime keeps the resources busy. Counts busy workers, not raw task spans, so the
# async PaRSEC GPU [submit,complete] spans (which overlap on a CUDA stream) don't
# push the curve above the worker count. Available for BOTH runtimes.
#   -> plots/sched_concurrency.{png,pdf}
#
# Usage:  plot_scheduler_health.r [base_dir]
#   default: data/manual_traces_*/runs. Works on any job with application.parquet
#   (StarPU) or tasks.parquet (PaRSEC).

suppressMessages({
  library(arrow); library(dplyr); library(tibble); library(ggplot2)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (nrow(runs) == 0) stop("no usable runs under ", base_dir)

# ---- (A) concurrency: busy WORKERS over time, app factorization kernels only ----
# Two corrections vs a naive task-span sweep:
#   - keep only FACTORIZATION_KERNELS (drop matrix-gen/init + sinks), the panel_st
#     compute tasks, so the curve is application work, not data generation.
#   - count busy WORKERS, not overlapping task spans: PaRSEC times a GPU task as
#     [submit, complete] (async), so several spans overlap on one CUDA stream and
#     a raw sweep exceeds the worker count. Unioning each worker's intervals first
#     caps every resource at one, keeping the curve within #workers.
concurrency <- function(run_dir) {
  ex <- read_exec(run_dir)                       # unified, normalised to ms
  if (is.null(ex) || nrow(ex) == 0) return(NULL)
  workers <- n_distinct(ex$worker)               # all execution resources (CPU+GPU)
  ex <- ex %>% filter(.data$kernel %in% FACTORIZATION_KERNELS)
  if (nrow(ex) == 0) return(NULL)
  busy <- ex %>%
    group_by(.data$worker) %>% arrange(.data$Start, .by_group = TRUE) %>%
    mutate(seg = cumsum(.data$Start > lag(cummax(.data$End), default = -Inf))) %>%
    group_by(.data$worker, .data$seg) %>%
    summarise(Start = min(.data$Start), End = max(.data$End), .groups = "drop")
  gmin <- min(busy$Start)
  bind_rows(
    tibble(t = busy$Start - gmin, d =  1L),
    tibble(t = busy$End   - gmin, d = -1L)
  ) %>%
    group_by(.data$t) %>% summarise(d = sum(.data$d), .groups = "drop") %>%
    arrange(.data$t) %>% mutate(running = cumsum(.data$d)) %>%
    mutate(workers = workers)
}

conc <- bind_rows(lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  c <- concurrency(r$path)
  if (is.null(c)) return(NULL)
  c$cfg <- paste0(r$runtime, ":", r$scheduler); c$algo <- r$algo; c
}))

if (!is.null(conc) && nrow(conc) > 0) {
  conc <- conc %>% mutate(panel = paste0(.data$algo, " ", .data$cfg))
  wk <- conc %>% group_by(.data$panel) %>% summarise(workers = first(.data$workers))
  p_conc <- ggplot(conc, aes(.data$t, .data$running)) +
    geom_step(aes(colour = .data$cfg), direction = "hv", linewidth = 0.4) +
    geom_hline(data = wk, aes(yintercept = .data$workers),
               linetype = "dashed", colour = "grey40") +
    facet_wrap(~panel, scales = "free_x") +
    labs(title = "Perfil de concorrencia: workers ocupados com computacao da aplicacao",
         subtitle = paste0(base_dir,
                          "  |  linha tracejada = no total de workers (CPU+GPU)"),
         x = "tempo (ms, t0 = 1o start)", y = "workers ocupados (tarefas de computacao)",
         colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  save_plot(p_conc, "sched_concurrency", width = 12, height = 7)
} else message("sem dados de concorrencia (nenhum application/tasks utilizavel)")
