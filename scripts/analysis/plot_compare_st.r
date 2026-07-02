#!/usr/bin/env Rscript
#
# Side-by-side StarVZ space-time (panel_st) for StarPU vs PaRSEC on the SAME
# problem, so the two runtimes' execution can be eyeballed and validated next to
# each other. This is the payoff of the parsec-cuda branch: the GPU PaRSEC runs
# now carry the full StarVZ phase-1 parquet set (application/dag/colors/y/...,
# reconstructed offline by parsec_phase1.sh), so the SAME native panel StarVZ
# draws for StarPU draws for PaRSEC, and they line up in one figure.
#
# For each algorithm found under base_dir it renders one composite figure laid
# out as a grid with StarPU on the LEFT column and PaRSEC on the RIGHT column,
# one scheduler per row:
#
#     starpu:dmda  | parsec:lfq
#     starpu:dmdas | parsec:gd
#
# Kernel-colour alignment: StarPU ("sgemm"/"spotrf"/"strsm"...) and PaRSEC
# ("gemm"/"potrf"/"Trsm"...) emit different Value labels with different
# colours.parquet hues, so the same op paints a different colour per runtime on
# a side-by-side panel. norm_st_kernels() (in trace_common.r) rewrites each
# run's Application$Value via norm_kernel() AND rebuilds Colors from the shared
# KERNEL_COLORS, so "gemm" lands on the same hue in both columns.
#
# Outputs <PLOTS_DIR>/st_compare_<algo>.{png,pdf}. PLOTS_DIR defaults to "plots"
# (honoured via trace_common.r); compare_all.sh sets it to plots/compare.
#
# Usage:  plot_compare_st.r [base_dir]
#   base_dir defaults to the GPU job (spotrf/sgeqrf n19200 b960, 4 StarPU + 4
#   PaRSEC runs). Needs each run's application.parquet (run parsec_phase1.sh for
#   the PaRSEC GPU runs first).

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

# panel_st knobs: same as plot_gpu_st.r -- aggregation on, ABE off (heterogeneous
# CPU+GPU set), idleness off (async GPU spans), CPB only when the run has a Dag.
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
  svz <- drop_init_kernels(svz)   # plgsy etc.: geracao da matriz, nao algoritmo
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

  # Build the by-row plot order so wrap_plots(ncol = 2) puts StarPU left,
  # PaRSEC right. Pad the shorter side with NULL so rows stay aligned.
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
