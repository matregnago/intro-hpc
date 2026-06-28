#!/usr/bin/env Rscript
#
# Extra panels (add-2, add-3 in the plan) to differentiate the schedulers:
#
#   (A) Concurrency profile -- number of compute tasks running simultaneously
#       over time, from application.parquet. Shows ramp-up / steady-state / tail
#       and how well each runtime keeps the cores busy. Available for BOTH
#       runtimes today (only needs application.parquet) -> matched comparison.
#       -> plots/sched_concurrency.{png,pdf}
#
#   (B) Ready / lack-ready (starvation) -- the number of eligible-but-not-running
#       tasks over time. When it hits 0 with free workers the scheduler is
#       starved. Read from variable.parquet (Type=="Ready") for StarPU; for PaRSEC
#       it is reconstructed from tasks.parquet (DAG + timing) via
#       parsec_ready_to_parquet.r. Runs without either source are skipped.
#       -> plots/sched_ready.{png,pdf}
#
# Usage:  plot_scheduler_health.r [base_dir]
#   default: the dpotrf n19200 b480 job (concurrency for all 4 runs; Ready for the
#   2 StarPU runs, whose variable.parquet carries it). Point at a job with PaRSEC
#   tasks.parquet (e.g. data/manual_traces_*/runs) for PaRSEC Ready curves.

suppressMessages({
  library(arrow); library(dplyr); library(tibble); library(ggplot2)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
# reconstruct_ready() for the PaRSEC path (skips its own CLI when sourced)
source(file.path(script_dir, "parsec_ready_to_parquet.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (nrow(runs) == 0) stop("no usable runs under ", base_dir)

# ---- (A) concurrency: sweep line over compute task [Start,End) (ms) ----
concurrency <- function(run_dir) {
  ex <- read_exec(run_dir)                       # unified, normalised to ms
  if (is.null(ex) || nrow(ex) == 0) return(NULL)
  gmin <- min(ex$Start)
  bind_rows(
    tibble(t = ex$Start - gmin, d =  1L),
    tibble(t = ex$End   - gmin, d = -1L)
  ) %>%
    group_by(.data$t) %>% summarise(d = sum(.data$d), .groups = "drop") %>%
    arrange(.data$t) %>% mutate(running = cumsum(.data$d)) %>%
    mutate(workers = n_distinct(ex$worker))
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
    labs(title = "Perfil de concorrencia: tarefas em execucao ao longo do tempo",
         subtitle = paste0(base_dir,
                          "  |  linha tracejada = no de workers (CPU)"),
         x = "tempo (ms, t0 = 1o start)", y = "tarefas em execucao",
         colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  save_plot(p_conc, "sched_concurrency", width = 12, height = 7)
} else message("sem dados de concorrencia (nenhum application/tasks utilizavel)")

# ---- (B) Ready / lack-ready over time ----
# StarPU: native variable.parquet Type=="Ready"; PaRSEC: reconstruct from tasks.
ready_series <- function(row) {
  vf <- file.path(row$path, "variable.parquet")
  if (row$variable && file.exists(vf)) {
    v <- read_parquet(vf) %>% filter(grepl("Ready", .data$Type))
    if (nrow(v) > 0) {
      g <- min(v$Start)
      return(tibble(Start = v$Start - g, Value = v$Value, src = "nativo"))
    }
  }
  if (row$tasks) {
    tk <- read_tasks_norm(row$path)              # ms
    v <- tryCatch(reconstruct_ready(tk), error = function(e) NULL)
    if (!is.null(v) && nrow(v) > 0) {
      return(tibble(Start = v$Start, Value = v$Value, src = "reconstruido"))
    }
  }
  NULL
}

ready <- bind_rows(lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  s <- ready_series(r)
  if (is.null(s)) { message("  sem Ready para ", r$dir); return(NULL) }
  s$cfg <- paste0(r$algo, " ", r$runtime, ":", r$scheduler); s
}))

if (!is.null(ready) && nrow(ready) > 0) {
  # starvation summary: share of timeline with Ready == 0
  starv <- ready %>% group_by(.data$cfg) %>%
    arrange(.data$Start, .by_group = TRUE) %>%
    mutate(dur = lead(.data$Start) - .data$Start) %>%
    filter(!is.na(.data$dur)) %>%
    summarise(span = sum(.data$dur),
              starved = sum(.data$dur[.data$Value <= 0]),
              mean_ready = weighted.mean(.data$Value, .data$dur),
              .groups = "drop") %>%
    mutate(starved_pct = 100 * .data$starved / .data$span)
  print(starv, n = Inf)
  dir.create("plots", showWarnings = FALSE)
  write.csv(starv, "plots/sched_ready_summary.csv", row.names = FALSE)

  p_ready <- ggplot(ready, aes(.data$Start, .data$Value, colour = .data$cfg)) +
    geom_step(direction = "hv", linewidth = 0.4) +
    facet_wrap(~cfg, scales = "free_x") +
    labs(title = "Tarefas prontas (Ready) ao longo do tempo -- starvation quando = 0",
         subtitle = paste0(base_dir,
                          "  |  Ready nativo (StarPU) ou reconstruido (PaRSEC)"),
         x = "tempo (us, t0 = 1o start)", y = "tarefas prontas (Ready)",
         colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  save_plot(p_ready, "sched_ready", width = 11, height = 6)
} else message("sem dados de Ready (nem variable nativo nem tasks p/ reconstruir)")
