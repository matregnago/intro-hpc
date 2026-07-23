#!/usr/bin/env Rscript
#
# Reconstruct StarPU's "Ready" scheduler counter for a PaRSEC run and write it
# as a StarVZ-style variable.parquet. PaRSEC's .prof never records the counter,
# but eligibility is reconstructable offline from tasks.parquet: a task is
# "ready" during [max(EndTime of its predecessors), StartTime), and
# Ready(t) is the step function (+1 at each ready_time, -1 at each StartTime).
#
# This is data-availability readiness (by the DAG), not the runtime queue's
# exact occupancy. The all-independent init tasks spike the counter at t0;
# --compute-only drops them for a factorization-only profile.
#
# Usage:  parsec_ready_to_parquet.r [--compute-only] <run_dir> [run_dir ...]
#   needs <run_dir>/tasks.parquet; writes <run_dir>/variable.parquet.

suppressMessages({
  library(arrow)
  library(dplyr)
  library(tibble)
})

# Kernels that are data init / sinks rather than factorization work.
INIT_SINK_KERNELS <- c("splgsy", "splrnt", "parsec_dtd_data_flush")

#' Reconstruct the Ready step-function from a tasks data frame.
reconstruct_ready <- function(tasks, compute_only = FALSE) {
  tk <- tasks %>% filter(!is.na(.data$StartTime), !is.na(.data$EndTime))
  if (compute_only) tk <- tk %>% filter(!(.data$Name %in% INIT_SINK_KERNELS))
  if (nrow(tk) == 0) stop("no timed tasks to reconstruct Ready from")

  endmap <- setNames(tk$EndTime, tk$JobId)
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

  bind_rows(
    tibble(t = ready_time,   d =  1L),
    tibble(t = tk$StartTime, d = -1L)
  ) %>%
    mutate(t = .data$t - gmin) %>%
    group_by(.data$t) %>%
    summarise(d = sum(.data$d), .groups = "drop") %>%
    arrange(.data$t) %>%
    mutate(Value = as.double(cumsum(.data$d))) %>%
    transmute(Start = .data$t, End = lead(.data$t), Value = .data$Value) %>%
    filter(!is.na(.data$End), .data$End > .data$Start) %>%
    mutate(
      Duration     = .data$End - .data$Start,
      Type         = "Ready",
      ResourceId   = "scheduler",     # panel_ready filters on "sched" + "Ready"
      Node         = 0L,
      Resource     = "scheduler",
      ResourceType = "scheduler"
    ) %>%
    select("ResourceId", "Type", "Start", "End", "Duration",
           "Value", "Node", "Resource", "ResourceType")
}

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

# CLI (skipped when this file is source()d).
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  compute_only <- "--compute-only" %in% args
  run_dirs <- args[args != "--compute-only"]
  if (length(run_dirs) == 0) {
    stop("usage: parsec_ready_to_parquet.r [--compute-only] <run_dir> [run_dir ...]")
  }
  for (d in run_dirs) write_ready_variable(d, compute_only = compute_only)
}
