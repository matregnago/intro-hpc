#!/usr/bin/env Rscript
#
# Side-by-side StarVZ space-time (panel_st) traces for StarPU vs PaRSEC on the
# SAME problem (dpotrf n19200 b480), so the two runtimes' execution can be
# compared visually. Both runs carry the full StarVZ phase-1 parquet set
# (application.parquet, ...); PaRSEC's came from dbp2paje + `starvz -1`.
#
# Outputs plots/trace_compare_st.{png,pdf} (StarPU on top, PaRSEC below) plus a
# standalone panel per run. Edit base_dir/runs when pointing at a new job.

suppressMessages({
  library(tidyverse); library(starvz); library(patchwork)
})

base_dir <- "data/parcial/chameleon-runtimes_782088/runs"
runs <- c(
  "StarPU (lws)" = "0001_starpu_lws_dpotrf_n19200_b480_rep1",
  "PaRSEC (lfq)" = "0008_parsec_lfq_dpotrf_n19200_b480_rep1"
)

configure <- function(svz) {
  svz$config$st$outliers          <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$step  <- 100
  svz$config$st$idleness          <- TRUE
  svz
}

slug <- function(s) str_replace_all(s, "[ ()/]+", "_")

panels <- list()
for (label in names(runs)) {
  path <- file.path(base_dir, runs[[label]])
  message("lendo: ", path)
  svz <- configure(starvz_read(path, selective = FALSE))
  svz <- norm_st_kernels(svz)
  p <- panel_st(svz) + ggtitle(label)
  panels[[label]] <- p
  ggsave(paste0("plots/trace_st_", slug(label), ".png"), p,
         width = 14, height = 5, dpi = 300, limitsize = FALSE)
}

dir.create("plots", showWarnings = FALSE)
combined <- panels[[1]] / panels[[2]] +
  plot_annotation(title = "Trace space-time: StarPU vs PaRSEC  (dpotrf n=19200, b=480)")
ggsave("plots/trace_compare_st.png", combined,
       width = 14, height = 10, dpi = 300, limitsize = FALSE)
ggsave("plots/trace_compare_st.pdf", combined,
       width = 14, height = 10, device = cairo_pdf, limitsize = FALSE)
message("salvo: plots/trace_compare_st.{png,pdf}")
