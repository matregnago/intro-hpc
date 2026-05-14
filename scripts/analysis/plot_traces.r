suppressMessages(library(tidyverse))
suppressMessages(library(starvz))
suppressMessages(library(patchwork))

base_dir <- "data/chameleon-runtimes_782088/runs"

runs <- c(
  "starpu/lws"  = "0001_starpu_lws_dpotrf_n19200_b480_rep1",
  "starpu/ws"   = "0006_starpu_ws_dpotrf_n19200_b480_rep1",
  "parsec/gd"   = "0007_parsec_gd_dpotrf_n19200_b480_rep1",
  "parsec/lfq"  = "0008_parsec_lfq_dpotrf_n19200_b480_rep1"
)

configure <- function(svz) {
  svz$config$st$outliers <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$step <- 100
  svz$config$st$idleness <- TRUE
  svz
}

slug <- function(label) gsub("/", "_", label, fixed = TRUE)

for (label in names(runs)) {
  dir_name <- runs[[label]]
  path <- file.path(base_dir, dir_name)
  message("lendo: ", path)

  svz <- starvz_read(path, selective = FALSE)
  svz <- configure(svz)
  plot <- panel_st(svz) + ggtitle(label)

  png_out <- paste0("plots/trace_panel_st_", slug(label), ".png")
  pdf_out <- paste0("plots/trace_panel_st_", slug(label), ".pdf")
  ggsave(png_out, plot,
    width = 16, height = 9, dpi = 300, limitsize = FALSE
  )
  ggsave(pdf_out, plot,
    width = 16, height = 9, device = cairo_pdf
  )
  message("salvo: ", png_out, " e ", pdf_out)
}
