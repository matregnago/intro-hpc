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
# Usage:  plot_kiteration.r <run_dir|base_dir> [run_dir ...]
#   Each arg is either a run dir (holding tasks.parquet) or a base dir (a job's
#   runs/ folder) that gets expanded to its run subdirs -- so this takes the SAME
#   single base_dir compare_all.sh hands every other comparison script. Needs
#   <run_dir>/tasks.parquet with a non-NA Iteration column (re-run
#   parsec_tasks_to_parquet.r if it predates the Iteration derivation).
#
#   Groups the runs BY ALGORITHM and, per algorithm, arranges every
#   scheduler/runtime in a 2-column grid on a SHARED time axis (within one algo
#   all runs share N/b, so makespans are comparable and the per-scheduler
#   "signatures" line up):
#     plots/kiteration_compare_<algo>.{png,pdf}   (>1 run for that algo)
#     plots/kiteration_<runtime>_<sched>_<algo>.{png,pdf}  (a lone run)
#   Row order is fixed (starpu:dmda, starpu:dmdas, parsec:lfq, parsec:gd) to match
#   plot_compare_st.r. Honours PLOTS_DIR (compare_all.sh routes here).

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
  message(sprintf("    %s [%s]: %d tasks over %d iterations",
                  label, meta$algo, nrow(app), max(app$Iteration) + 1L))
  list(panel = p, label = label,
       runtime   = ifelse(is.na(meta$runtime), "run", meta$runtime),
       scheduler = ifelse(is.na(meta$scheduler), basename(run_dir), meta$scheduler),
       algo      = ifelse(is.na(meta$algo), "run", meta$algo),
       x_max = max(app$End), niter = max(app$Iteration) + 1L)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: plot_kiteration.r <run_dir|base_dir> [run_dir ...]")

# Each arg is a run dir (has tasks.parquet) or a base dir to expand to its runs.
expand_arg <- function(a) {
  if (file.exists(file.path(a, "tasks.parquet"))) return(a)
  subs <- list.dirs(a, recursive = FALSE)
  subs[file.exists(file.path(subs, "tasks.parquet"))]
}
run_dirs <- unique(unlist(lapply(args, expand_arg)))
if (length(run_dirs) == 0)
  stop("no run dirs with tasks.parquet under: ", paste(args, collapse = ", "))

panels <- Filter(Negate(is.null), lapply(run_dirs, build_panel))
if (length(panels) == 0) stop("no runs with usable k-iteration data")

# Fixed scheduler/runtime row order, matching plot_compare_st.r's grid.
CONFIG_ORDER <- c("starpu:dmda", "starpu:dmdas", "parsec:lfq", "parsec:gd")
config_rank <- function(p) {
  m <- match(paste0(p$runtime, ":", p$scheduler), CONFIG_ORDER)
  if (is.na(m)) length(CONFIG_ORDER) + 1L else m
}
stem_of <- function(p) paste0("kiteration_", p$runtime, "_", p$scheduler, "_", p$algo)

# One figure per algorithm: stack all schedulers/runtimes on a SHARED time axis
# (within an algo every run shares N/b, so makespans are comparable) and collect
# the (identical) kernel legend once. A lone run for an algo -> single panel.
by_algo <- split(panels, vapply(panels, function(p) p$algo, character(1)))
for (algo in names(by_algo)) {
  grp <- by_algo[[algo]]
  grp <- grp[order(vapply(grp, config_rank, integer(1)))]
  if (length(grp) == 1L) {
    pn <- grp[[1]]
    save_plot(pn$panel + ggtitle(paste0("k-iteration  ", pn$label, "  ", algo)),
              stem_of(pn), width = 14, height = 8)
    next
  }
  x_end <- max(vapply(grp, function(p) p$x_max, numeric(1)))
  # 2-column landscape grid (fits a 16:9 slide), same config order as
  # plot_compare_st.r; the shared xlim keeps makespans comparable.
  combined <- wrap_plots(
    lapply(grp, function(p) p$panel + coord_cartesian(xlim = c(0, x_end))),
    ncol = 2
  ) +
    plot_layout(guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom"))
  save_plot(combined, paste0("kiteration_compare_", algo),
            width = 16, height = 4.5 * ceiling(length(grp) / 2))
}
