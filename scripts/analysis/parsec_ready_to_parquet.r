#!/usr/bin/env Rscript
#
# Reconstruct StarPU's "Ready" scheduler counter for a PaRSEC run and write it as
# a StarVZ-style variable.parquet -- the same table StarVZ builds from StarPU's
# Paje "Variable" events (PajeSetVariable). PaRSEC's .prof never records this
# counter (the scheduling-event traces in parsec/scheduling.c are commented out),
# so dbp2paje yields an empty variable.parquet. But "Ready" is the number of
# tasks that are eligible-but-not-yet-running, which is RECONSTRUCTABLE offline
# from the DAG + per-task timing we already extract -- no runtime rebuild needed.
#
# How:
#   - input is the run's tasks.parquet (from parsec_tasks_to_parquet.r): each row
#     is a task with StartTime/EndTime and DependsOn = space-separated predecessor
#     JobIds (multi-valued; in a Cholesky a gemm has 2-3 preds, trsm/syrk 1-2).
#   - a task T becomes ready at  ready_time(T) = max(EndTime of its predecessors)
#     (entry tasks, with no recorded preds, become ready at the global start), and
#     leaves the ready set when it starts executing (StartTime). So T is "ready"
#     during [ready_time(T), StartTime(T)).
#   - Ready(t) = #{T : ready_time(T) <= t < StartTime(T)}, built as a step function
#     (+1 at each ready_time, -1 at each StartTime, cumulative sum over time).
#
# Output schema = StarVZ variable.parquet:
#   ResourceId, Type, Start, End, Duration, Value, Node, Resource, ResourceType
# with ResourceId="scheduler" (panel_ready filters grepl("sched", ResourceId) &
# grepl("Ready", Type)). Times are shifted so the earliest start is 0.
#
# Caveats:
#   - this is data-availability readiness (eligibility by the DAG), not the
#     runtime scheduler's exact queue occupancy (it ignores scheduling latency).
#   - the init tasks (splgsy/splrnt) are all independent, so they spike the
#     counter at t0 (a real backlog). Pass --compute-only to drop the init/sink
#     kernels and get the factorization-only profile.
#   - StarPU's "Submitted" counter is NOT reconstructable this way (it needs the
#     real DTD insertion timestamps, which are not in the trace).
#
# Usage:
#   parsec_ready_to_parquet.r [--compute-only] <run_dir> [run_dir ...]
# Each run_dir must contain tasks.parquet; writes <run_dir>/variable.parquet.

suppressMessages({
  library(arrow)
  library(dplyr)
  library(tibble)
})

# Kernels that are data init / sinks rather than factorization work.
INIT_SINK_KERNELS <- c("splgsy", "splrnt", "parsec_dtd_data_flush")

#' Reconstruct the Ready step-function from a tasks data frame.
#' @param tasks       data frame with JobId, Name, StartTime, EndTime, DependsOn.
#' @param compute_only drop INIT_SINK_KERNELS before reconstructing.
#' @return a StarVZ-schema variable.parquet tibble (Type = "Ready").
reconstruct_ready <- function(tasks, compute_only = FALSE) {
  tk <- tasks %>% filter(!is.na(.data$StartTime), !is.na(.data$EndTime))
  if (compute_only) tk <- tk %>% filter(!(.data$Name %in% INIT_SINK_KERNELS))
  if (nrow(tk) == 0) stop("no timed tasks to reconstruct Ready from")

  endmap <- setNames(tk$EndTime, tk$JobId)   # JobId -> EndTime
  gmin   <- min(tk$StartTime)

  # ready_time(T) = max(EndTime of timed predecessors); entry tasks -> gmin.
  ready_time <- vapply(seq_len(nrow(tk)), function(i) {
    dep <- tk$DependsOn[i]
    if (is.na(dep) || !nzchar(dep)) {
      return(gmin)
    }
    pe <- endmap[strsplit(dep, " ", fixed = TRUE)[[1]]]
    pe <- pe[!is.na(pe)]
    if (!length(pe)) gmin else max(pe)
  }, numeric(1))
  ready_time <- pmin(ready_time, tk$StartTime)   # causality clamp

  # step function: +1 when a task becomes ready, -1 when it starts running.
  bind_rows(
    tibble(t = ready_time,   d =  1L),
    tibble(t = tk$StartTime, d = -1L)
  ) %>%
    mutate(t = .data$t - gmin) %>%               # shift so t0 = first start
    group_by(.data$t) %>%
    summarise(d = sum(.data$d), .groups = "drop") %>%
    arrange(.data$t) %>%
    mutate(Value = as.double(cumsum(.data$d))) %>%
    transmute(Start = .data$t, End = lead(.data$t), Value = .data$Value) %>%
    filter(!is.na(.data$End), .data$End > .data$Start) %>%
    mutate(
      Duration     = .data$End - .data$Start,
      Type         = "Ready",
      ResourceId   = "scheduler",
      Node         = 0L,
      Resource     = "scheduler",
      ResourceType = "scheduler"
    ) %>%
    select("ResourceId", "Type", "Start", "End", "Duration",
           "Value", "Node", "Resource", "ResourceType")
}

#' Read a run's tasks.parquet, reconstruct Ready, write variable.parquet into it.
write_ready_variable <- function(run_dir, compute_only = FALSE) {
  tasks_file <- file.path(run_dir, "tasks.parquet")
  if (!file.exists(tasks_file)) {
    warning("no tasks.parquet in ", run_dir, " -- skipping")
    return(invisible(NULL))
  }
  v <- reconstruct_ready(read_parquet(tasks_file), compute_only = compute_only)
  out <- file.path(run_dir, "variable.parquet")
  write_parquet(v, out)
  message(sprintf(
    "%s: Ready max=%g mean=%.1f (%d intervals) -> %s",
    run_dir, max(v$Value), weighted.mean(v$Value, v$Duration), nrow(v), out
  ))
  invisible(v)
}

# CLI (skipped when this file is source()d, e.g. by plot_ready.r).
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  compute_only <- "--compute-only" %in% args
  run_dirs <- args[args != "--compute-only"]
  if (length(run_dirs) == 0) {
    stop("usage: parsec_ready_to_parquet.r [--compute-only] <run_dir> [run_dir ...]")
  }
  for (d in run_dirs) write_ready_variable(d, compute_only = compute_only)
}
