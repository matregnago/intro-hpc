#!/usr/bin/env Rscript
#
# Shared helpers for the runtime-comparison plots: parse run-dir names, discover
# runs, normalise kernel names across StarPU/PaRSEC, and save figures.
# Sourced, never run.

suppressMessages({
  library(dplyr); library(stringr); library(tibble); library(arrow)
})

# Stable hue per normalised kernel stem, so the same op keeps the same colour
# across runs/runtimes/schedulers.
KERNEL_COLORS <- c(
  potrf = "#e41a1c", trsm = "#377eb8", syrk = "#4daf4a", gemm = "#984ea3",
  herk = "#4daf4a", trmm = "#377eb8",
  geqrt = "#e41a1c", ormqr = "#984ea3", qr_couple = "#377eb8", qr_apply = "#4daf4a",
  getrf = "#e41a1c", getrf_nopiv = "#e41a1c"
)

#' StarVZ-style Colors table for a set of kernel Values: stable hue via
#' KERNEL_COLORS (looked up by the aliased stem), deterministic hue_pal
#' fallback for unknown kernels.
KERNEL_COLORS_TO_TIBBLE <- function(values) {
  vals  <- sort(unique(as.character(values)))
  stems <- norm_kernel(vals, raw = FALSE)
  cols  <- unname(KERNEL_COLORS[stems])
  miss  <- is.na(cols)
  if (any(miss)) cols[miss] <- scales::hue_pal()(sum(miss))
  tibble(Value = factor(vals, levels = vals), Color = cols, Use = TRUE)
}

#' Drop matrix-generation/init tasks from a starvz_read() result's Application.
#' StarPU FxT traces don't contain them, so on a side-by-side ST panel they only
#' show up on the PaRSEC column. Call BEFORE norm_st_kernels().
drop_init_kernels <- function(svz) {
  if (is.null(svz$Application) || !"Value" %in% names(svz$Application)) return(svz)
  keep <- !norm_kernel(svz$Application$Value) %in% INIT_KERNELS
  svz$Application <- svz$Application[keep, ]
  svz
}

#' Normalise Application$Value and rebuild Colors of a starvz_read() result, so
#' panel_st paints the same kernel the same hue for StarPU ("dgemm") and PaRSEC
#' ("gemm").
norm_st_kernels <- function(svz) {
  if (is.null(svz$Application) || !"Value" %in% names(svz$Application)) return(svz)
  app <- svz$Application
  app$Value <- norm_kernel(app$Value)
  if ("lowercase" %in% names(app))
    app$lowercase <- tolower(as.character(app$Value))
  svz$Application <- app
  svz$Colors <- KERNEL_COLORS_TO_TIBBLE(app$Value)
  svz
}

# Compute-kernel stems (precision-stripped). NB the Chameleon QR apply codelets
# end in 't' (tsmqrt/ttmqrt/tpmqrt): StarPU emits the TS/TT scheme, PaRSEC the
# unified TP scheme; KERNEL_ALIAS lines the two up.
CHOLESKY_KERNELS <- c("potrf", "trsm", "syrk", "gemm", "herk", "trmm")
QR_KERNELS       <- c("geqrt", "ormqr", "unmqr",
                      "tsqrt", "ttqrt", "tpqrt",        # couple factor
                      "tsmqrt", "ttmqrt", "tpmqrt")     # couple apply
LU_KERNELS       <- c("getrf", "getrf_nopiv")
INIT_KERNELS     <- c("plgsy", "plrnt", "plghe", "lacpy", "laset")
COMPUTE_KERNELS  <- c(CHOLESKY_KERNELS, QR_KERNELS, LU_KERNELS, INIT_KERNELS,
                      "qr_couple", "qr_apply")

# Same math under different codelet names: ormqr = unmqr; StarPU's TS/TT and
# PaRSEC's TP tile-QR schemes collapse to one label per QR stage.
KERNEL_ALIAS <- c(
  unmqr  = "ormqr",
  tsqrt  = "qr_couple",  ttqrt  = "qr_couple",  tpqrt  = "qr_couple",
  tsmqrt = "qr_apply",   ttmqrt = "qr_apply",   tpmqrt = "qr_apply"
)

#' Collapse a kernel state name to a common label: strip the d/s/c/z precision
#' prefix ("dgemm"->"gemm"), then apply KERNEL_ALIAS. raw=TRUE skips the alias
#' step (e.g. to tell tsqrt from tpqrt in a per-runtime breakdown).
norm_kernel <- function(value, raw = FALSE) {
  v <- tolower(as.character(value))
  stem <- ifelse(
    nchar(v) > 1L &
      substr(v, 1, 1) %in% c("d", "s", "c", "z") &
      substr(v, 2, nchar(v)) %in% COMPUTE_KERNELS,
    substr(v, 2, nchar(v)), v
  )
  if (raw) return(stem)
  ifelse(stem %in% names(KERNEL_ALIAS), KERNEL_ALIAS[stem], stem)
}

#' Parse a run-dir basename "<idx>_<runtime>_<sched>_<algo>_n<N>_b<B>_rep<R>"
#' (the scheme scripts/run.sh writes) into its experiment factors; fields that
#' don't match come back NA.
parse_run_id <- function(dir_name) {
  d <- basename(dir_name)
  m <- str_match(
    d,
    "^(\\d+)_([a-z]+)_([a-z]+)_([a-z][a-z0-9_]*)_n(\\d+)_b(\\d+)_rep(\\d+)"
  )
  tibble(
    dir       = d,
    idx       = m[, 2],
    runtime   = m[, 3],
    scheduler = m[, 4],
    algo      = m[, 5],
    n         = as.integer(m[, 6]),
    b         = as.integer(m[, 7]),
    rep       = as.integer(m[, 8]),
    label     = ifelse(is.na(m[, 3]), d,
                       paste0(m[, 3], ":", m[, 4]))
  )
}

#' Which of the comparison-relevant parquets a run dir holds.
run_has <- function(run_dir) {
  has <- function(f) file.exists(file.path(run_dir, f))
  list(
    application = has("application.parquet"),
    tasks       = has("tasks.parquet"),
    dag         = has("dag.parquet"),
    variable    = has("variable.parquet"),
    states      = has("states.parquet")
  )
}

#' Discover runs under a base dir (a job's runs/ directory): one row per run
#' with parse_run_id() columns + `path` + the run_has() availability flags.
list_runs <- function(base_dir, pattern = NULL) {
  dirs <- list.dirs(base_dir, recursive = FALSE)
  if (!is.null(pattern)) dirs <- dirs[str_detect(basename(dirs), pattern)]
  if (length(dirs) == 0) {
    stop("no run dirs under ", base_dir,
         if (!is.null(pattern)) paste0(" matching /", pattern, "/") else "")
  }
  meta <- bind_rows(lapply(dirs, parse_run_id))
  meta$path <- dirs
  avail <- bind_rows(lapply(dirs, function(d) as_tibble(run_has(d))))
  bind_cols(meta, avail) %>% arrange(.data$dir)
}

#' Scale that brings a run's dag/tasks/states times to MILLISECONDS. The offline
#' PaRSEC converters -- the only writers of a states.parquet -- emit
#' microseconds; StarVZ phase-1 parquets are already ms.
time_scale_to_ms <- function(run_dir) {
  if (file.exists(file.path(run_dir, "states.parquet"))) 1e-3 else 1
}

#' Read a run's tasks.parquet with StartTime/EndTime normalised to ms.
read_tasks_norm <- function(run_dir) {
  f <- file.path(run_dir, "tasks.parquet")
  if (!file.exists(f)) return(NULL)
  s <- time_scale_to_ms(run_dir)
  read_parquet(f) %>%
    mutate(StartTime = .data$StartTime * s, EndTime = .data$EndTime * s)
}

#' Unified per-kernel execution-interval table for a run, runtime-agnostic:
#' prefer application.parquet, fall back to tasks.parquet (the source that is
#' GPU-complete for PaRSEC). Times in ms.
#' @return tibble(kernel, Start, End, Duration, worker, src) of COMPUTE_KERNELS
#'   only, or NULL if neither source has usable compute states.
read_exec <- function(run_dir) {
  pick <- function(tb) tb %>%
    filter(.data$kernel %in% COMPUTE_KERNELS, !is.na(.data$Duration),
           .data$Duration > 0)

  appf <- file.path(run_dir, "application.parquet")
  if (file.exists(appf)) {
    a <- read_parquet(appf)
    if (nrow(a) > 0 && "Value" %in% names(a)) {
      out <- pick(tibble(kernel = norm_kernel(a$Value), Start = a$Start,
                         End = a$End, Duration = a$Duration,
                         worker = as.character(a$ResourceId), src = "application"))
      if (nrow(out) > 0) return(out)
    }
  }
  tkf <- file.path(run_dir, "tasks.parquet")
  if (file.exists(tkf)) {
    s <- time_scale_to_ms(run_dir)
    t <- read_parquet(tkf)
    out <- pick(tibble(kernel = norm_kernel(t$Name), Start = t$StartTime * s,
                       End = t$EndTime * s, Duration = (t$EndTime - t$StartTime) * s,
                       worker = as.character(t$WorkerId), src = "tasks"))
    if (nrow(out) > 0) return(out)
  }
  NULL
}

#' Output directory for figures/tables; honours $PLOTS_DIR (default "plots").
plots_dir <- function() {
  d <- Sys.getenv("PLOTS_DIR", "plots")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

#' Save PNG (300 dpi) + PDF (cairo) under plots_dir().
save_plot <- function(plot, stem, width = 12, height = 6) {
  od  <- plots_dir()
  png <- file.path(od, paste0(stem, ".png"))
  pdf <- file.path(od, paste0(stem, ".pdf"))
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 300,
                  limitsize = FALSE)
  ggplot2::ggsave(pdf, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, limitsize = FALSE)
  message("salvo: ", png, " e ", pdf)
}
