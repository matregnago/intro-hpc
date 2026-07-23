#!/usr/bin/env Rscript
#
# Side-by-side StarVZ space-time (panel_st) for StarPU vs PaRSEC on the same
# problem. One composite figure per algorithm: StarPU on the left column,
# PaRSEC on the right, one scheduler per row. norm_st_kernels() aligns the
# kernel names/colours across the two runtimes.
#
# Usage:  plot_compare_st.r [base_dir]
#   needs each run's application.parquet (run parsec_phase1.sh for the PaRSEC
#   GPU runs first). Outputs <PLOTS_DIR>/st_compare_<algo>.{png,pdf}.

suppressMessages({
  library(tidyverse); library(starvz); library(patchwork)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"

# panel_st knobs: aggregation on, ABE off (heterogeneous CPU+GPU set),
# idleness off (async GPU spans), CPB only when the run has a Dag.
configure <- function(svz) {
  svz$config$st$outliers           <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$step   <- 100
  svz$config$st$idleness           <- FALSE
  svz$config$st$makespan           <- TRUE
  svz$config$st$abe$active         <- FALSE
  svz$config$st$abe$label          <- FALSE
  svz$config$st$tasks$active       <- FALSE
  svz$config$st$cpb <- !is.null(svz$Dag) && nrow(svz$Dag) > 0
  svz
}

panel_for <- function(run_dir, label) {
  message("lendo: ", run_dir)
  svz <- configure(starvz_read(run_dir, selective = FALSE))
  svz <- drop_init_kernels(svz)
  svz <- norm_st_kernels(svz)
  panel_st(svz) + ggtitle(label)
}

runs <- list_runs(base_dir) %>% filter(.data$application)
if (nrow(runs) == 0) stop("no runs with application.parquet under ", base_dir,
                          " (run parsec_phase1.sh for the PaRSEC runs)")

for (alg in sort(unique(runs$algo))) {
  sub <- runs %>% filter(.data$algo == alg) %>% arrange(.data$runtime, .data$scheduler)
  starpu <- sub %>% filter(.data$runtime == "starpu")
  parsec <- sub %>% filter(.data$runtime == "parsec")
  if (nrow(starpu) == 0 || nrow(parsec) == 0) {
    message("skip ", alg, ": need both runtimes (have starpu=", nrow(starpu),
            ", parsec=", nrow(parsec), ")")
    next
  }

  # Row-major plot order so wrap_plots(ncol = 2) puts StarPU left, PaRSEC
  # right; the shorter side is padded to keep rows aligned.
  nrows <- max(nrow(starpu), nrow(parsec))
  cell  <- function(df, i) if (i <= nrow(df))
    panel_for(df$path[i], paste0(df$runtime[i], ":", df$scheduler[i]))
  else plot_spacer()
  plots <- list()
  for (i in seq_len(nrows)) {
    plots[[length(plots) + 1]] <- cell(starpu, i)
    plots[[length(plots) + 1]] <- cell(parsec, i)
  }

  n  <- unique(na.omit(sub$n)); b <- unique(na.omit(sub$b))
  title <- sprintf("Space-time: StarPU vs PaRSEC  -  %s (n=%s, b=%s)",
                   alg, paste(n, collapse = "/"), paste(b, collapse = "/"))
  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(title = title)

  save_plot(combined, paste0("st_compare_", alg),
            width = 22, height = 5.2 * nrows)
}
