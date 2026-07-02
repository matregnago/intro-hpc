#!/usr/bin/env Rscript
#
# Render the stock StarVZ space-time panel (panel_st) for a PaRSEC GPU run, using
# the application/y/colors/dag/variable parquets reconstructed offline by
# scripts/analysis/parsec_phase1.sh. This is the payoff of the synthetic
# Application slot: the SAME native panel StarVZ draws for StarPU now draws for
# PaRSEC on GPU, where the official dbp2paje phase 1 produces no application.parquet.
#
# Overlays enabled:
#   - CPB (critical-path bound) when the run carries a dag.parquet (it does) --
#     hl_global_cpb() needs only $Dag + $Application sharing JobId.
#   - makespan.
# Overlays kept OFF:
#   - ABE (geom_abe): heterogeneous CPU+GPU resource set -> the homogeneous area
#     bound is invalid; the valid GPU floor is in plot_bounds.r.
#   - idleness (geom_idleness): PaRSEC times a GPU task as [submit, complete]
#     (async), so spans on one CUDA stream overlap and per-lane "idle" goes
#     negative/undefined. Meaningful only on the CPU lanes, so left OFF to avoid
#     misleading labels on the CUDA lanes.
#   - task-dependency arrows (st$tasks): need the StarPU $Last table we don't have.
#
# Optionally also emits panel_imbalance and panel_utilheatmap (same Application,
# no extra data) when --extra is passed.
#
# Usage:
#   plot_gpu_st.r [--extra] [run_dir ...]
#   Defaults to the 4 PaRSEC GPU runs under data/manual_traces_20260621_172820.
#   Writes plots/gpu_st_<runtime>_<sched>_<algo>.{png,pdf} (one per run), plus
#   plots/gpu_{imbalance,utilheatmap}_<...> with --extra.

suppressMessages({
  library(tidyverse); library(starvz)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))   # parse_run_id, save_plot

args  <- commandArgs(trailingOnly = TRUE)
extra <- "--extra" %in% args
runs  <- args[args != "--extra"]
if (length(runs) == 0) {
  runs <- Sys.glob("data/manual_traces_20260621_172820/runs/*_parsec_*")
}
if (length(runs) == 0) stop("no run dirs given and no default PaRSEC GPU runs found")

# Same knobs as plot_traces.r's configure(), but ABE off (GPU) and CPB gated on Dag.
configure <- function(svz) {
  svz$config$st$outliers           <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$step   <- 100
  svz$config$st$idleness           <- FALSE   # async GPU spans -> idle undefined
  svz$config$st$makespan           <- TRUE
  svz$config$st$abe$active         <- FALSE
  svz$config$st$abe$label          <- FALSE
  svz$config$st$tasks$active       <- FALSE
  svz$config$st$cpb <- !is.null(svz$Dag) && nrow(svz$Dag) > 0
  svz
}

for (rd in runs) {
  rd <- sub("/$", "", rd)
  if (!file.exists(file.path(rd, "application.parquet"))) {
    message("skip (no application.parquet -- run parsec_phase1.sh first): ", rd)
    next
  }
  meta  <- parse_run_id(rd)
  label <- if (is.na(meta$runtime)) basename(rd) else
    paste0(meta$runtime, ":", meta$scheduler, "  ", meta$algo)
  stem  <- if (is.na(meta$runtime)) paste0("gpu_st_", basename(rd)) else
    paste0("gpu_st_", meta$runtime, "_", meta$scheduler, "_", meta$algo)
  message("lendo: ", rd)

  svz <- configure(starvz_read(rd, selective = FALSE))
  svz <- drop_init_kernels(svz)   # plgsy etc.: geracao da matriz, nao algoritmo
  svz <- norm_st_kernels(svz)
  p   <- panel_st(svz) + ggtitle(label)
  save_plot(p, stem, width = 16, height = 9)

  if (extra) {
    tag <- sub("^gpu_st_", "", stem)
    save_plot(panel_imbalance(svz)  + ggtitle(paste0("imbalance  ", label)),
              paste0("gpu_imbalance_", tag),  width = 14, height = 4)
    save_plot(panel_utilheatmap(svz) + ggtitle(paste0("util heatmap  ", label)),
              paste0("gpu_utilheatmap_", tag), width = 14, height = 5)
  }
}
