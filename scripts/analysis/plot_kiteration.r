#!/usr/bin/env Rscript
#
# k-iteration plot (iterations on Y, time on X) via StarVZ's panel_kiteration(),
# to compare the schedulers' "signatures" across the outer factorization loop.
# PaRSEC has no Application$Iteration in its trace; parsec_tasks_to_parquet.r
# derives it, and panel_kiteration only needs $Application/$Colors/$config, so
# a hand-built starvz_data list is enough (no starvz_read).
#
# Usage:  plot_kiteration.r <run_dir|base_dir> [run_dir ...]
#   run dirs need tasks.parquet with a non-NA Iteration column.
# Output: one figure per algorithm, all schedulers/runtimes on a shared time
#   axis: plots/kiteration_compare_<algo>.{png,pdf} (or kiteration_<cfg>_<algo>
#   for a lone run). Honours PLOTS_DIR.

suppressMessages({
  library(ggplot2); library(starvz); library(patchwork)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

# Colors table in panel_kiteration's expected shape (Value character, no Use),
# with the shared KERNEL_COLORS hues + hue_pal fallback.
build_colors <- function(values) {
  vals <- sort(unique(values))
  cols <- KERNEL_COLORS[vals]
  missing <- is.na(cols)
  if (any(missing)) cols[missing] <- scales::hue_pal()(sum(missing))
  tibble(Value = vals, Color = unname(cols))
}

# panel_kiteration for one run; NULL if the run has no usable data.
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
      Value        = norm_kernel(.data$Name)
    ) %>%
    # factorization work only: names(KERNEL_COLORS) is exactly the post-alias
    # factorization kernel set (matrix-gen and runtime sinks drop out).
    filter(.data$Value %in% names(KERNEL_COLORS))
  if (nrow(app) == 0) { message("no timed iterated tasks in ", run_dir, ", skip"); return(NULL) }

  svz <- structure(list(
    Application = app,
    Colors      = build_colors(app$Value),
    config = list(
      base_size = 22, expand = 0.05,
      limits = list(start = NA, end = NA),
      # subite/pernode must be explicit FALSE: panel_kiteration errors on NULL.
      kiteration = list(subite = FALSE, pernode = FALSE,
                        middlelines = NULL, legend = TRUE)
    )
  ), class = "starvz_data")

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

# Fixed row order, matching plot_compare_st.r's grid.
CONFIG_ORDER <- c("starpu:dmda", "starpu:dmdas", "parsec:lfq", "parsec:gd")
config_rank <- function(p) {
  m <- match(paste0(p$runtime, ":", p$scheduler), CONFIG_ORDER)
  if (is.na(m)) length(CONFIG_ORDER) + 1L else m
}
stem_of <- function(p) paste0("kiteration_", p$runtime, "_", p$scheduler, "_", p$algo)

# One figure per algorithm: within an algo every run shares N/b, so a shared
# time axis keeps makespans comparable.
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
  combined <- wrap_plots(
    lapply(grp, function(p) p$panel + coord_cartesian(xlim = c(0, x_end))),
    ncol = 2
  ) +
    plot_layout(guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom"))
  save_plot(combined, paste0("kiteration_compare_", algo),
            width = 16, height = 4.5 * ceiling(length(grp) / 2))
}
