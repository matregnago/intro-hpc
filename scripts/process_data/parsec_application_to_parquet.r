#!/usr/bin/env Rscript
#
# Build the StarVZ tables the official PaRSEC phase 1 can't produce on a GPU
# run: application.parquet + y.parquet + colors.parquet. With them present,
# starvz_read() fills $Application/$Y/$Colors and the stock StarVZ panels run
# on PaRSEC GPU traces.
#
# Source is the run's tasks.parquet (from parsec_tasks_to_parquet.r): one row
# per timed task, with a JobId that joins dag.parquet and a ResourceId
# ("CPU<N>" / "CUDA<dev>_<stream>") that decides its lane and ResourceType.
# Times converted us -> ms to match StarPU's clock.
#
# Usage:  parsec_application_to_parquet.r <run_dir> [run_dir ...]
#   needs <run_dir>/tasks.parquet. Normally invoked by parsec_phase1.sh.

suppressMessages({
  library(arrow); library(dplyr); library(stringr); library(tibble)
})

# norm_kernel() + KERNEL_COLORS_TO_TIBBLE() for the colors table.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "..", "analysis", "trace_common.r"))

# Column order (+ types) of StarVZ's StarPU application.parquet (25 cols).
APP_COLS <- c(
  "ResourceId", "Start", "End", "Duration", "Depth", "Value", "Size", "Params",
  "Footprint", "Tag", "JobId", "SubmitOrder", "GFlop", "X", "Y", "Iteration",
  "Subiteration", "NumaNodes", "Node", "Resource", "ResourceType", "lowercase",
  "Outlier", "Height", "Position"
)

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

  # Y: one lane per resource, mirroring StarVZ's StarPU layout so both runtimes
  # read the same on a side-by-side panel_st: CUDA streams first (bottom of the
  # panel, double height), then CPU workers in descending index order.
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
      Start    = .data$StartTime / 1000,                   # us -> ms
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

  colors <- KERNEL_COLORS_TO_TIBBLE(vals)

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
