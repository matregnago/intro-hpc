#!/usr/bin/env Rscript
#
# Idea #3 (graficos.md): the k-iteration plot -- iterations on Y, time on X --
# to see whether each scheduler leaves a different "signature" across the outer
# loop of a right-looking factorization.
#
# StarVZ ships exactly this as panel_kiteration() (R/phase2_kchart.R), but it
# needs an Application$Iteration column that StarPU's trace carries and PaRSEC's
# does NOT. We recover k for PaRSEC in parsec_tasks_to_parquet.r (one diagonal
# kernel per iteration over the DTD submit order) and feed StarVZ's own panel a
# hand-built starvz_data list -- panel_kiteration only touches $Application,
# $Colors and a few $config fields, so no application.parquet / starvz_read is
# needed. This is the GPU-complete path (tasks.parquet, not dbp2paje).
#
# Usage:  plot_kiteration.r <run_dir> [run_dir ...]
#   needs <run_dir>/tasks.parquet with a non-NA Iteration column (re-run
#   parsec_tasks_to_parquet.r if it predates the Iteration derivation).
#   Writes plots/kiteration_<runtime>_<sched>.{png,pdf} per run, and -- when given
#   more than one run -- a stacked comparison plots/kiteration_compare.{png,pdf}
#   with a SHARED time axis (so the per-scheduler "signatures" line up).

suppressMessages({
  library(ggplot2); library(starvz); library(patchwork)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

# Stable colors for the usual factorization kernels, so a kernel keeps its hue
# across runs/runtimes; anything else gets filled in from the default hue palette.
KERNEL_COLORS <- c(
  potrf = "#e41a1c", trsm = "#377eb8", syrk = "#4daf4a", gemm = "#984ea3",
  herk = "#4daf4a", trmm = "#377eb8",
  geqrt = "#e41a1c", ormqr = "#984ea3", qr_couple = "#377eb8", qr_apply = "#4daf4a"
)

build_colors <- function(values) {
  vals <- sort(unique(values))
  cols <- KERNEL_COLORS[vals]
  missing <- is.na(cols)
  if (any(missing)) cols[missing] <- scales::hue_pal()(sum(missing))
  tibble(Value = vals, Color = unname(cols))
}

# Build the k-iteration panel for one run via StarVZ's panel_kiteration(). Returns
# a list(panel, label, stem, x_max, niter) or NULL if the run has no usable data.
build_panel <- function(run_dir) {
  meta <- parse_run_id(run_dir)
  tasks <- read_tasks_norm(run_dir)            # StartTime/EndTime -> ms
  if (is.null(tasks)) { message("no tasks.parquet in ", run_dir, ", skip"); return(NULL) }
  if (!"Iteration" %in% names(tasks) || all(is.na(tasks$Iteration))) {
    message("no usable Iteration in ", run_dir,
            " (re-run parsec_tasks_to_parquet.r), skip"); return(NULL)
  }

  app <- tasks %>%
    filter(!is.na(.data$StartTime), !is.na(.data$EndTime), !is.na(.data$Iteration)) %>%
    transmute(
      Iteration    = as.integer(.data$Iteration),
      Subiteration = 0L,
      Node         = 0L,
      Start        = .data$StartTime,
      End          = .data$EndTime,
      Value        = norm_kernel(.data$Name)   # precision-agnostic stem
    ) %>%
    # keep only factorization work: matrix-gen (splgsy) lands before the first
    # opener (Iteration NA, already dropped); runtime data-flush tasks would
    # otherwise show up against the last k. names(KERNEL_COLORS) is exactly the
    # post-alias factorization kernel set.
    filter(.data$Value %in% names(KERNEL_COLORS))
  if (nrow(app) == 0) { message("no timed iterated tasks in ", run_dir, ", skip"); return(NULL) }

  svz <- structure(list(
    Application = app,
    Colors      = build_colors(app$Value),
    config = list(
      base_size = 22, expand = 0.05,
      limits = list(start = NA, end = NA),
      # subite=FALSE -> Y is the real Iteration; pernode=FALSE -> no Node facet.
      # Both must be explicit: panel_kiteration does `if (per_node)` (errors on NULL).
      kiteration = list(subite = FALSE, pernode = FALSE,
                        middlelines = NULL, legend = TRUE)
    )
  ), class = "starvz_data")   # panel_kiteration's starvz_check_data requires this

  label <- meta$label
  p <- panel_kiteration(svz) + ggtitle(label) + ylab("Iteration (k)")
  stem <- paste0("kiteration_",
                 ifelse(is.na(meta$runtime), "run", meta$runtime), "_",
                 ifelse(is.na(meta$scheduler), basename(run_dir), meta$scheduler))
  message(sprintf("    %s: %d tasks over %d iterations",
                  label, nrow(app), max(app$Iteration) + 1L))
  list(panel = p, label = label, stem = stem,
       x_max = max(app$End), niter = max(app$Iteration) + 1L,
       runtime = ifelse(is.na(meta$runtime), "run", meta$runtime))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: plot_kiteration.r <run_dir> [run_dir ...]")

panels <- Filter(Negate(is.null), lapply(args, build_panel))
if (length(panels) == 0) stop("no runs with usable k-iteration data")

# Per-run figures (own time scale).
for (pn in panels) save_plot(pn$panel + ggtitle(paste0("k-iteration  ", pn$label)),
                             pn$stem, width = 14, height = 8)

# Comparison figure: stack all panels on a SHARED time axis so the per-scheduler
# signatures are directly comparable; collect the (identical) kernel legend once.
if (length(panels) > 1) {
  x_end <- max(vapply(panels, function(p) p$x_max, numeric(1)))
  combined <- wrap_plots(
    lapply(panels, function(p) p$panel + coord_cartesian(xlim = c(0, x_end))),
    ncol = 1
  ) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "k-iteration por escalonador (eixo de tempo compartilhado)",
                    theme = theme(legend.position = "bottom"))
  rt <- panels[[1]]$runtime
  save_plot(combined, paste0("kiteration_compare_", rt),
            width = 14, height = 5 * length(panels))
}
