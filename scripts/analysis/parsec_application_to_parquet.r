#!/usr/bin/env Rscript
#
# Build the StarVZ slot that the official PaRSEC phase 1 (dbp2paje + `starvz -1`)
# CANNOT produce on a GPU run: application.parquet (+ its companions y.parquet and
# colors.parquet). With these three present, starvz_read() fills $Application/$Y/
# $Colors and the stock StarVZ panels that key off the Application slot run on
# PaRSEC GPU traces -- panel_st (space-time), the CPB critical-path overlay
# (needs $Dag + $Application sharing JobId, both of which we already have), and
# panel_imbalance / panel_utilheatmap / panel_progress / panel_resource_usage_task.
#
# Source is the run's tasks.parquet (from parsec_tasks_to_parquet.r): one row per
# DAG task with a real JobId that JOINS the dag.parquet, the WorkerId it ran on,
# StartTime/EndTime (microseconds) and the derived Iteration. We prefer it over
# states.parquet because states unions overlapping intervals and loses per-task
# identity (so its JobId can't be matched against the critical path).
#
# What each output table is, and why:
#   application.parquet -- one row per timed task, in StarVZ's 25-column StarPU
#       Application schema. Times converted us->ms to match StarPU's clock. This is
#       a HYBRID run: tasks carry a ResourceId from the tasks converter ("CPU<N>"
#       for the 16 CPU workers, "CUDA<dev>_<stream>" for the GPU streams), so we
#       build one lane per resource and set ResourceType CPU/CUDA accordingly --
#       the GPU-executed kernels (often most of the gemms) land on their own CUDA
#       lanes instead of being dropped. ABE is still left OFF in plot_gpu_st.r (the
#       homogeneous area bound is invalid on a heterogeneous CPU+GPU set; the valid
#       GPU floor lives in plot_bounds.r).
#   y.parquet           -- one row per ResourceId (Parent="CPU<N>"/"CUDA<d>_<s>",
#       Type="Worker State", Height=1, Position=0..n-1), the StarPU Y format;
#       yconf() turns Parent into the axis ResourceId. CPUs first, then CUDA streams.
#   colors.parquet      -- one Color per distinct kernel Value (stable hues from
#       KERNEL_COLORS via the precision-agnostic stem, hue_pal fallback).
#
# Usage:  parsec_application_to_parquet.r <run_dir> [run_dir ...]
#   needs <run_dir>/tasks.parquet. Writes application/y/colors.parquet there.
#   Normally invoked by scripts/analysis/parsec_phase1.sh after the tasks/dag
#   converters; runnable standalone too.

suppressMessages({
  library(arrow); library(dplyr); library(stringr); library(tibble)
})

# norm_kernel() (precision-agnostic stem + QR aliasing) for the color lookup.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

# Stable hues for the usual factorization kernels (mirrors plot_kiteration.r), so a
# kernel keeps its color across runs/runtimes; looked up by the aliased stem so
# spotrf/dpotrf->potrf, stpqrt->qr_couple, etc. all land on the same hue.
KERNEL_COLORS <- c(
  potrf = "#e41a1c", trsm = "#377eb8", syrk = "#4daf4a", gemm = "#984ea3",
  herk = "#4daf4a", trmm = "#377eb8",
  geqrt = "#e41a1c", ormqr = "#984ea3", qr_couple = "#377eb8", qr_apply = "#4daf4a"
)

# Final column order (+ types) of StarVZ's StarPU application.parquet (25 cols).
APP_COLS <- c(
  "ResourceId", "Start", "End", "Duration", "Depth", "Value", "Size", "Params",
  "Footprint", "Tag", "JobId", "SubmitOrder", "GFlop", "X", "Y", "Iteration",
  "Subiteration", "NumaNodes", "Node", "Resource", "ResourceType", "lowercase",
  "Outlier", "Height", "Position"
)

#' Colors table for a set of (raw) kernel Values: stable hue by aliased stem,
#' hue_pal fallback for anything unknown. Value is a factor with the given levels.
build_colors <- function(values) {
  vals  <- sort(unique(as.character(values)))
  stems <- norm_kernel(vals, raw = FALSE)        # spotrf->potrf, stpqrt->qr_couple
  cols  <- unname(KERNEL_COLORS[stems])
  miss  <- is.na(cols)
  if (any(miss)) cols[miss] <- scales::hue_pal()(sum(miss))
  tibble(Value = factor(vals, levels = vals), Color = cols, Use = TRUE)
}

convert_one <- function(run_dir) {
  tf <- file.path(run_dir, "tasks.parquet")
  if (!file.exists(tf)) {
    message("  no tasks.parquet in ", run_dir, " (run parsec_tasks_to_parquet.r first), skip")
    return(invisible())
  }
  message("==> ", run_dir)

  tk <- read_parquet(tf) %>%
    filter(!is.na(.data$StartTime), !is.na(.data$EndTime), !is.na(.data$ResourceId))
  if (nrow(tk) == 0) {
    message("  tasks.parquet has no timed tasks with a ResourceId, skip")
    return(invisible())
  }

  # Y: one lane per resource present, StarPU "Worker State" format, mirroring the
  # layout StarVZ phase 1 emits for StarPU so both runtimes read the same on a
  # side-by-side panel_st: CUDA streams FIRST (Position 0, 2, 4, ... -- the
  # bottom of the panel) with Height=2 (double-height GPU lanes), highest stream
  # index first; then the CPU workers with Height=1 in DESCENDING index order
  # (CPU<max> right above the CUDA block, CPU0 at the top).
  res      <- unique(tk$ResourceId)
  is_cpu   <- grepl("^CPU", res)
  cpu_ord  <- res[is_cpu][order(as.integer(sub("^CPU", "", res[is_cpu])),
                                decreasing = TRUE)]
  cuda_ord <- sort(res[!is_cpu], decreasing = TRUE)
  lanes    <- c(cuda_ord, cpu_ord)
  heights  <- c(rep(2, length(cuda_ord)), rep(1, length(cpu_ord)))
  y <- tibble(
    Parent   = lanes,
    Type     = "Worker State",
    Height   = heights,
    Position = cumsum(c(0, heights[-length(heights)]))
  )
  posmap <- setNames(y$Position, y$Parent)
  hgtmap <- setNames(y$Height,   y$Parent)

  vals <- sort(unique(as.character(tk$Name)))

  app <- tk %>%
    mutate(
      Start    = .data$StartTime / 1000,                   # us -> ms (StarPU clock)
      End      = .data$EndTime   / 1000,
      Duration = (.data$EndTime - .data$StartTime) / 1000
    ) %>%
    transmute(
      ResourceId   = .data$ResourceId,
      Start        = .data$Start,
      End          = .data$End,
      Duration     = .data$Duration,
      Depth        = 0,
      Value        = factor(as.character(.data$Name), levels = vals),
      Size         = 0L,
      Params       = NA_character_,
      Footprint    = NA_character_,
      Tag          = NA_character_,
      JobId        = as.character(.data$JobId),            # joins dag.parquet
      SubmitOrder  = as.character(.data$SubmitOrder),
      GFlop        = NA_real_,
      X            = NA_integer_,
      Y            = NA_integer_,
      Iteration    = as.integer(.data$Iteration),
      Subiteration = 0L,
      NumaNodes    = NA_character_,
      Node         = factor("0"),
      Resource     = factor(.data$ResourceId),
      ResourceType = factor(if_else(grepl("^CUDA", .data$ResourceId), "CUDA", "CPU")),
      lowercase    = tolower(as.character(.data$Name)),
      Outlier      = FALSE,
      Height       = unname(hgtmap[.data$ResourceId]),
      Position     = unname(posmap[.data$ResourceId])
    ) %>%
    select(all_of(APP_COLS))

  colors <- build_colors(vals)

  write_parquet(app,    file.path(run_dir, "application.parquet"))
  write_parquet(y,      file.path(run_dir, "y.parquet"))
  write_parquet(colors, file.path(run_dir, "colors.parquet"))

  message(sprintf(
    "    application=%d rows  lanes=%d (CPU=%d, CUDA=%d)  kernels={%s}  span=%.1fms -> application/y/colors.parquet",
    nrow(app), length(lanes), sum(is_cpu), sum(!is_cpu),
    paste(vals, collapse = ","), max(app$End) - min(app$Start)
  ))
  invisible(app)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: parsec_application_to_parquet.r <run_dir> [run_dir ...]")
for (d in args) convert_one(d)
