#!/usr/bin/env Rscript
#
# Shared helpers for the StarVZ-based runtime-comparison plots (plot_task_times.r,
# plot_bounds.r, plot_wait_time.r, plot_scheduler_health.r). Centralises the three
# things every comparison needs and that differ between StarPU and PaRSEC:
#
#   1. parse_run_id()  -- turn a run-dir name like
#        "0008_parsec_lfq_dpotrf_n19200_b480_rep1" into (idx,runtime,scheduler,
#        algo,n,b,rep) so runs can be faceted/labelled uniformly.
#   2. norm_kernel()   -- StarPU names kernels dgemm/dsyrk/dtrsm/dpotrf and PaRSEC
#        gemm/syrk/Trsm/potrf; both collapse to a precision-agnostic stem
#        (gemm/syrk/trsm/potrf/...) so the same kernel lines up across runtimes.
#   3. list_runs()     -- discover the run dirs under a job's runs/ directory and
#        attach the parsed metadata + which StarVZ parquets each one actually has
#        (the PaRSEC dag/tasks tables exist only where the grapher .dot was
#        collected and reconstructed -- see parsec_*_to_parquet.r).
#
# Sourced, not run. Pure helpers; no side effects.

suppressMessages({
  library(dplyr); library(stringr); library(tibble); library(arrow)
})

# Factorization + init/util kernels we treat as "compute" task types. Anything
# not here (StarPU's B/Callback/Fi/P, PaRSEC's Waiting/movein/moveout/cuda and the
# scheduler/runtime classes) is runtime bookkeeping, dropped from task-time stats.
CHOLESKY_KERNELS <- c("potrf", "trsm", "syrk", "gemm", "herk", "trmm")
# Tile-QR kernels. NB the Chameleon codelet names end in 't' for the apply
# kernels (tsmqrt/ttmqrt/tpmqrt), which an earlier version misspelled (tsmqr...)
# and silently dropped -- they are the bulk of QR work. StarPU emits the classic
# TS/TT scheme (tsqrt+ttqrt / tsmqrt+ttmqrt); PaRSEC emits the unified TP scheme
# (tpqrt / tpmqrt). See KERNEL_ALIAS for lining the two up.
QR_KERNELS       <- c("geqrt", "ormqr", "unmqr",
                      "tsqrt", "ttqrt", "tpqrt",        # couple factor
                      "tsmqrt", "ttmqrt", "tpmqrt")     # couple apply
INIT_KERNELS     <- c("plgsy", "plrnt", "plghe", "lacpy", "laset")  # data init / copy
# Valid kernels = the raw stems (used by norm_kernel's prefix-strip check) PLUS
# the aliased labels norm_kernel can emit, so filtering on its (aliased) output
# keeps the QR couple/apply work instead of dropping it.
COMPUTE_KERNELS  <- c(CHOLESKY_KERNELS, QR_KERNELS, INIT_KERNELS,
                      "qr_couple", "qr_apply")

# Map runtime-specific kernel names onto a common operation, so the SAME math
# lines up across runtimes in the comparisons. `ormqr`/`unmqr` are literally the
# same op under two codelet names. The TS/TT (StarPU) and TP (PaRSEC) tile-QR
# schemes are different decompositions of the same two steps -- grouped to one
# logical label each so both runtimes have the same QR stages (the per-task mean
# then mixes StarPU's two sub-kernels; the totals/counts stay comparable).
KERNEL_ALIAS <- c(
  unmqr  = "ormqr",
  tsqrt  = "qr_couple",  ttqrt  = "qr_couple",  tpqrt  = "qr_couple",
  tsmqrt = "qr_apply",   ttmqrt = "qr_apply",   tpmqrt = "qr_apply"
)

#' Collapse a runtime-specific kernel state name to a common operation label.
#' Strips the d/s/c/z precision prefix ("dgemm"->"gemm", "stsmqrt"->"tsmqrt"),
#' then applies KERNEL_ALIAS so equivalent ops (ormqr=unmqr, TS/TT=TP) share a
#' label. Unknown names come back lower-cased and unchanged.
#' @param raw if TRUE, skip KERNEL_ALIAS and keep the precision-stripped kernel
#'   name (e.g. to tell tsqrt from tpqrt in a per-runtime breakdown).
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

#' Parse a run-dir basename into its experiment factors.
#' Expects "<idx>_<runtime>_<sched>_<algo>_n<N>_b<B>_rep<R>" (the scheme
#' scripts/run.sh writes). Fields that don't match come back NA.
parse_run_id <- function(dir_name) {
  d <- basename(dir_name)
  m <- str_match(
    d,
    "^(\\d+)_([a-z]+)_([a-z]+)_([a-z]+)_n(\\d+)_b(\\d+)_rep(\\d+)"
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
    # compact, human label for facets/titles
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

#' Discover runs under a base dir (a job's runs/ directory) and attach metadata.
#' @param base_dir directory containing the per-run subdirectories.
#' @param pattern  optional regex to keep only some run dirs.
#' @return a tibble: one row per run with parse_run_id() columns + `path` + the
#'   run_has() availability flags (application/tasks/dag/variable/states).
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

#' Time-unit scale to bring a run's tables to MILLISECONDS.
#' StarVZ-produced parquets (StarPU FxT, or PaRSEC via dbp2paje) are already in ms
#' and always carry an application.parquet. The PaRSEC adapter
#' (parsec_*_to_parquet.r) emits microseconds and is the ONLY source when there is
#' no application.parquet. So: application present -> ms (x1); tasks-only -> us
#' (x1e-3). This keeps StarPU and PaRSEC on the same clock when compared.
time_scale_to_ms <- function(run_dir) {
  if (file.exists(file.path(run_dir, "application.parquet"))) 1 else 1e-3
}

#' Total GPU device-busy time for a run, in MILLISECONDS, or NA if the run has no
#' GPU activity (so this doubles as the GPU-vs-CPU discriminator). This is the
#' valid "area" bound on the GPU, where the homogeneous ABE (Sum dur / nworkers)
#' breaks: PaRSEC's tasks.parquet pins all GPU work to ~16 fictional CPU thread
#' workers (no per-task device tag), so ABE can exceed the makespan. The device
#' actually executes its kernels back-to-back, so Sum(device-busy) is a single-GPU
#' throughput floor that IS below the makespan.
#'   - PaRSEC adapter: states.parquet rows typed "GPU" (the `cuda` device spans).
#'   - StarPU: application.parquet rows with ResourceType == "CUDA".
#' NB this is a conservative floor: with >1 concurrent CUDA stream the spans can
#' overlap, so the sum may overcount; we divide StarPU's by the #CUDA workers, and
#' for a single GPU it stays below the makespan on the validated traces.
gpu_busy <- function(run_dir) {
  stf <- file.path(run_dir, "states.parquet")
  if (file.exists(stf)) {
    st <- read_parquet(stf)
    if ("Type" %in% names(st)) {
      g <- st %>% filter(.data$Type == "GPU")
      if (nrow(g) > 0) return(sum(g$Duration) * time_scale_to_ms(run_dir))
    }
  }
  appf <- file.path(run_dir, "application.parquet")
  if (file.exists(appf)) {
    a <- read_parquet(appf)
    if (all(c("ResourceType", "Duration") %in% names(a))) {
      g <- a %>% filter(.data$ResourceType == "CUDA")
      if (nrow(g) > 0) return(sum(g$Duration) / n_distinct(g$ResourceId))
    }
  }
  NA_real_
}

#' Read a run's tasks.parquet with StartTime/EndTime normalised to ms.
read_tasks_norm <- function(run_dir) {
  f <- file.path(run_dir, "tasks.parquet")
  if (!file.exists(f)) return(NULL)
  s <- time_scale_to_ms(run_dir)
  read_parquet(f) %>%
    mutate(StartTime = .data$StartTime * s, EndTime = .data$EndTime * s)
}

#' Unified per-kernel execution-interval table for a run, runtime-agnostic.
#' StarPU and PaRSEC expose their work in different parquets, and for the GPU
#' traces the PaRSEC application.parquet can't be built (dbp2paje/pj_dump fails on
#' the GPU paje), while tasks.parquet (reconstructed from .dot + .prof) IS
#' GPU-complete. So we prefer application.parquet (both runtimes have it for the
#' CPU job) and fall back to tasks.parquet (the source that exists for PaRSEC on
#' the GPU job). Both measure the same thing -- a compute kernel's execution span.
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
    s <- time_scale_to_ms(run_dir)              # PaRSEC adapter us -> ms
    t <- read_parquet(tkf)
    out <- pick(tibble(kernel = norm_kernel(t$Name), Start = t$StartTime * s,
                       End = t$EndTime * s, Duration = (t$EndTime - t$StartTime) * s,
                       worker = as.character(t$WorkerId), src = "tasks"))
    if (nrow(out) > 0) return(out)
  }
  NULL
}

#' Standardised PNG+PDF save (300 dpi PNG + cairo_pdf), mirroring the existing
#' plot scripts. Creates plots/ if needed.
save_plot <- function(plot, stem, width = 12, height = 6) {
  dir.create("plots", showWarnings = FALSE)
  png <- paste0("plots/", stem, ".png")
  pdf <- paste0("plots/", stem, ".pdf")
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 300,
                  limitsize = FALSE)
  ggplot2::ggsave(pdf, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, limitsize = FALSE)
  message("salvo: ", png, " e ", pdf)
}
