#!/usr/bin/env Rscript
#
# Build a StarVZ-style tasks.parquet for a PaRSEC run from its grapher .dot +
# binary .prof, mirroring the per-task metadata table StarVZ derives from
# StarPU's tasks.rec. One row per task (DAG node), carrying its dependencies,
# the worker it ran on, and its timing.
#
# What maps from PaRSEC (and what doesn't):
#   - JobId/Name/Priority/{tpid}      <- .dot node: tooltip + label
#       label is "<thread/vp> name(tid)[tid]<priority>{tpid}"
#   - DependsOn / DepLabels           <- .dot incoming edges (src + edge label)
#   - StartTime / EndTime / WorkerId  <- .prof via dbp2xml (per-task span + the
#                                        THREAD it executed on); tpid:tid joins
#                                        structure to timing (see dag converter)
#   - SubmitOrder                     <- rank over (tpid, tid); DTD assigns task
#                                        ids at insertion, so this is submit order
#   - Iteration                       <- DERIVED, not in the trace. Chameleon submits
#                                        a right-looking factorization, so each outer
#                                        iteration k opens with exactly one diagonal
#                                        kernel (potrf/geqrt/getrf). Counting those
#                                        over SubmitOrder gives the exact k; tasks
#                                        before the first opener (matrix gen) -> NA.
#                                        This feeds StarVZ's panel_kiteration().
#   - Model/File/Line/Footprint/GFlop/Parameters/SubmitTime/Tag/
#     NumaNodes/Control               <- not in the PaRSEC trace -> NA
#       (GFlop is derivable analytically per kernel+tile; left NA for now.)
#
# Also writes <run_dir>/states.parquet: every timed interval the .prof records
# (not just compute kernels) as a StarVZ-style long table, one row per event,
# categorised into Computation / Scheduling / Idle / Runtime / Memory /
# DataTransfer / GPU. parsec.nix builds with PARSEC_PROF_TRACE_SCHEDULING_EVENTS
# and _ACTIVE_ARENA_SET, so the scheduler/idle/memory spans are already in the
# trace -- this surfaces them for the per-worker occupancy breakdown that
# distinguishes the schedulers (lfq vs gd). dbp2xml runs once for both tables.
#
# Usage:  parsec_tasks_to_parquet.r <run_dir> [run_dir ...]
#   needs cham_*.dot + cham_*.prof-* in run_dir, dbp2xml on PATH or $DBP2XML.
#   Writes <run_dir>/tasks.parquet and <run_dir>/states.parquet. See
#   parsec_dag_to_parquet.r for the companion dag.parquet and the shared
#   dbp2xml/out.xml gotcha.

suppressMessages({
  library(xml2); library(arrow); library(dplyr); library(stringr); library(tidyr)
})

# PaRSEC tags every timed interval with a dictionary "class". Compute kernels
# (gemm/potrf/...) become task rows; everything else is a runtime/scheduler/
# memory state the build emits (PARSEC_PROF_TRACE_SCHEDULING_EVENTS /
# _ACTIVE_ARENA_SET). We keep them all (categorised) in states.parquet for the
# per-worker occupancy breakdown. Anything not listed -- i.e. a kernel -- falls
# through to "Computation".
event_category <- function(name) {
  dplyr::case_when(
    name %in% c("Sched POLL", "Sched PUSH", "Queue REMOVE")          ~ "Scheduling",
    name == "Sched SLEEP"                                            ~ "Idle",
    name %in% c("RELEASE_DEPS", "ACTIVATE_CB", "DATA_FLUSH",
                "Device delegate", "PUT_CB")                         ~ "Runtime",
    name %in% c("MEMALLOC", "ARENA_MEMORY", "ARENA_ACTIVE_SET",
                "TASK_MEMORY")                                       ~ "Memory",
    name %in% c("movein", "moveout")                                 ~ "DataTransfer",
    name == "cuda"                                                   ~ "GPU",
    startsWith(name, "MPI_")                                         ~ "Communication",
    TRUE                                                             ~ "Computation"
  )
}

# Diagonal kernels that open one outer iteration k of a right-looking
# factorization (Cholesky/QR/LU). Exactly one fires per k, so counting them over
# the DTD submission order recovers the iteration index (see Iteration below).
ITER_OPENERS <- c("potrf", "geqrt", "getrf")

# Final schema (column order + types) of StarVZ's StarPU tasks.parquet.
TASKS_COLS <- c(
  "Control", "JobId", "SubmitOrder", "SubmitTime", "MPIRank", "Name", "Model",
  "File", "Line", "Priority", "DependsOn", "DepLabels", "Tag", "WorkerId",
  "MemoryNode", "StartTime", "EndTime", "Footprint", "GFlop", "Iteration",
  "Parameters", "NumaNodes"
)

# Long "state" table: one row per occupancy interval (overlapping recordings of
# the same class on a worker are unioned, see build_states), anchored to the
# same t0 as tasks so both tables share a clock. Times in microseconds.
STATES_COLS <- c("Worker", "ResourceId", "Type", "Value",
                 "Start", "End", "Duration", "JobId")

resolve_dbp2xml <- function() {
  env <- Sys.getenv("DBP2XML")
  if (nzchar(env) && file.exists(env)) return(env)
  found <- Sys.which("dbp2xml")
  if (nzchar(found)) return(unname(found))
  message("dbp2xml not on PATH; building .#parsec (first run only)...")
  out <- suppressWarnings(system2(
    "nix", c("build", "--impure", "--no-link", "--print-out-paths", ".#parsec"),
    stdout = TRUE, stderr = FALSE
  ))
  out <- tail(out[nzchar(out)], 1)
  bin <- if (length(out)) file.path(out, "bin", "dbp2xml") else ""
  if (nzchar(bin) && file.exists(bin)) return(bin)
  stop("dbp2xml not found: put it on PATH, set $DBP2XML, or run from a tree ",
       "where `nix build --impure .#parsec` works")
}

# .dot: per-node id/kernel/priority + edges (predecessor -> task, with label)
parse_dot <- function(dot_path) {
  ln <- readLines(dot_path, warn = FALSE)

  decl <- str_match(
    ln,
    paste0('^\\s*(\\S+)\\s+\\[.*<(-?\\d+)>\\{\\d+\\}".*',
           'tooltip="tpid=(\\d+):did=\\d+:tname=([^:"]+):tid=(\\d+)"')
  )
  decl <- decl[!is.na(decl[, 1]), , drop = FALSE]
  nodes <- tibble(
    name     = decl[, 2],
    Priority = as.integer(decl[, 3]),
    key      = paste(decl[, 4], decl[, 6], sep = ":"),   # tpid:tid
    Name     = decl[, 5]
  ) |> distinct(name, .keep_all = TRUE)
  name2key <- setNames(nodes$key, nodes$name)

  em <- str_match(ln, '^\\s*(\\S+)\\s*->\\s*(\\S+)\\s*\\[label="([^"]*)"')
  em <- em[!is.na(em[, 1]), , drop = FALSE]
  edges <- tibble(
    src   = unname(name2key[em[, 2]]),
    dst   = unname(name2key[em[, 3]]),
    label = em[, 4]
  ) |> filter(!is.na(src), !is.na(dst))

  list(nodes = nodes, edges = edges)
}

# .prof -> dbp2xml -> one row per timed event across ALL dictionary classes:
# the class name, the tpid:tid key it is attached to, the THREAD it ran on, and
# its [start,end] in ns. Both the per-task timing and the states table are
# derived from this, so dbp2xml runs only once.
parse_prof_events <- function(prof_path, dbp2xml) {
  prof_abs <- normalizePath(prof_path)
  work <- tempfile("dbp2xml_"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old <- setwd(work); on.exit(setwd(old), add = TRUE)
  # dbp2xml is a self-contained Nix binary (carries its own RPATH); the analysis
  # dev shell's R wrapper exports an LD_LIBRARY_PATH with a different glibc that
  # makes it fail to load, so run the child with LD_LIBRARY_PATH cleared.
  system2(dbp2xml, shQuote(prof_abs), stdout = FALSE, stderr = FALSE,
          env = "LD_LIBRARY_PATH=")
  xml_path <- file.path(work, "out.xml")
  if (!file.exists(xml_path) || file.info(xml_path)$size == 0) {
    stop("dbp2xml produced no out.xml for ", prof_path)
  }
  doc <- read_xml(xml_path)

  keys <- xml_find_all(doc, "//DICTIONARY/KEY")
  id2name <- setNames(xml_text(xml_find_first(keys, "./NAME")),
                      xml_attr(keys, "ID"))

  rows <- lapply(xml_find_all(doc, "//THREAD"), function(th) {
    ident  <- xml_text(xml_find_first(th, "./IDENTIFIER"))
    worker <- as.integer(str_match(ident, "Thread (\\d+)")[, 2])
    blocks <- xml_find_all(th, "./KEY")
    if (length(blocks) == 0) return(NULL)
    bind_rows(lapply(blocks, function(b) {
      evs <- xml_find_all(b, "./EVENT")
      if (length(evs) == 0) return(NULL)
      tibble(
        class    = unname(id2name[xml_attr(b, "ID")]),
        key      = xml_text(xml_find_first(evs, "./ID")),
        start_ns = as.numeric(xml_text(xml_find_first(evs, "./START"))),
        end_ns   = as.numeric(xml_text(xml_find_first(evs, "./END"))),
        WorkerId = worker
      )
    }))
  })
  rows <- bind_rows(rows)
  if (nrow(rows) == 0) {
    return(tibble(class = character(), key = character(), start_ns = numeric(),
                  end_ns = numeric(), WorkerId = integer()))
  }
  rows
}

# compute kernels only, collapsed to one [start,end] + worker per task key
# (min start / max end on the earliest worker -- matches the original timing).
compute_timing <- function(ev) {
  ev <- ev[event_category(ev$class) == "Computation", , drop = FALSE]
  if (nrow(ev) == 0) {
    return(tibble(key = character(), start_ns = numeric(),
                  end_ns = numeric(), WorkerId = integer()))
  }
  ev |>
    group_by(key) |>
    summarise(start_ns = min(start_ns), end_ns = max(end_ns),
              WorkerId = WorkerId[which.min(start_ns)], .groups = "drop")
}

# every timed interval as a long StarVZ-style state row, anchored to t0.
# PaRSEC records some events more than once with overlapping spans (each compute
# kernel appears twice), so we union overlapping intervals within each
# (worker, class): occupancy is not double-counted, while genuinely separate
# events -- gaps between them -- stay as distinct rows.
build_states <- function(ev, t0) {
  empty <- tibble(Worker = integer(), ResourceId = character(),
                  Type = character(), Value = character(), Start = numeric(),
                  End = numeric(), Duration = numeric(), JobId = character())
  ev <- ev[!is.na(ev$WorkerId), , drop = FALSE]
  if (nrow(ev) == 0) return(empty)
  out <- ev |>
    arrange(WorkerId, class, start_ns, end_ns) |>
    group_by(WorkerId, class) |>
    # new interval whenever this event starts after the running max end so far
    mutate(.grp = cumsum(start_ns > lag(cummax(end_ns), default = -Inf))) |>
    group_by(WorkerId, class, .grp) |>
    summarise(start_ns = min(start_ns), end_ns = max(end_ns),
              JobId = first(key), .groups = "drop") |>
    transmute(
      Worker     = as.integer(WorkerId),
      ResourceId = paste0("CPU", WorkerId),
      Type       = event_category(class),
      Value      = class,
      Start      = (start_ns - t0) / 1000,
      End        = (end_ns   - t0) / 1000,
      Duration   = (end_ns - start_ns) / 1000,
      JobId      = JobId
    ) |>
    arrange(Worker, Start) |>
    select(all_of(STATES_COLS))
  if (nrow(out) == 0) empty else out
}

convert_one <- function(run_dir, dbp2xml) {
  dot  <- list.files(run_dir, pattern = "\\.dot$", full.names = TRUE)
  prof <- list.files(run_dir, pattern = "\\.prof-", full.names = TRUE)
  prof <- prof[!grepl("\\.tar\\.gz$", prof)]
  if (length(dot) == 0)  { message("  no .dot in ", run_dir, ", skip");  return(invisible()) }
  if (length(prof) == 0) { message("  no .prof in ", run_dir, ", skip"); return(invisible()) }

  message("==> ", run_dir)
  g  <- parse_dot(dot[1])
  ev <- parse_prof_events(prof[1], dbp2xml)
  t  <- compute_timing(ev)

  t0 <- if (nrow(t)) min(t$start_ns) else 0
  timing <- t |> transmute(
    key, WorkerId,
    StartTime = (start_ns - t0) / 1000, EndTime = (end_ns - t0) / 1000
  )

  # collapse incoming edges into space-separated DependsOn / DepLabels
  deps <- g$edges |>
    arrange(dst, src) |>
    group_by(key = dst) |>
    summarise(DependsOn = paste(src, collapse = " "),
              DepLabels = paste(label, collapse = " "), .groups = "drop")

  tasks <- g$nodes |>
    transmute(key, Name, Priority,
              tpid = as.integer(str_extract(key, "^\\d+")),
              tid  = as.integer(str_extract(key, "\\d+$"))) |>
    arrange(tpid, tid) |>
    mutate(SubmitOrder = row_number()) |>
    # k = (#openers seen so far) - 1, over submission order. Tasks before the
    # first opener (matrix generation) get -1 -> NA, so they drop out cleanly.
    mutate(Iteration = {
      it <- cumsum(tolower(Name) %in% ITER_OPENERS) - 1L
      ifelse(it < 0L, NA_integer_, it)
    }) |>
    left_join(deps, by = "key") |>
    left_join(timing, by = "key") |>
    transmute(
      Control     = NA_character_,
      JobId       = key,
      SubmitOrder = as.integer(SubmitOrder),
      SubmitTime  = NA_real_,
      MPIRank     = 0L,
      Name        = Name,
      Model       = NA_character_,
      File        = NA_character_,
      Line        = NA_real_,
      Priority    = as.integer(Priority),
      DependsOn   = DependsOn,
      DepLabels   = DepLabels,
      Tag         = NA_character_,
      WorkerId    = as.integer(WorkerId),
      MemoryNode  = 0L,
      StartTime   = StartTime,
      EndTime     = EndTime,
      Footprint   = NA_character_,
      GFlop       = NA_real_,
      Iteration   = as.integer(Iteration),
      Parameters  = NA_character_,
      NumaNodes   = NA_character_
    ) |>
    select(all_of(TASKS_COLS))

  out <- file.path(run_dir, "tasks.parquet")
  write_parquet(tasks, out)

  n_timed <- sum(!is.na(tasks$StartTime))
  niter <- suppressWarnings(max(tasks$Iteration, na.rm = TRUE)) + 1L
  message(sprintf(
    "    tasks=%d (timed=%d, %.0f%%) iters=%s workers=%s prio=[%d..%d] -> %s",
    nrow(tasks), n_timed, 100 * n_timed / nrow(tasks),
    if (is.finite(niter)) niter else "NA",
    paste(sort(unique(na.omit(tasks$WorkerId))), collapse = ","),
    min(tasks$Priority, na.rm = TRUE), max(tasks$Priority, na.rm = TRUE), out
  ))

  states <- build_states(ev, t0)
  states_out <- file.path(run_dir, "states.parquet")
  write_parquet(states, states_out)

  brk <- states |>
    group_by(Type) |>
    summarise(ms = sum(Duration) / 1000, .groups = "drop") |>
    arrange(desc(ms))
  message(sprintf(
    "    states=%d -> %s\n      %s", nrow(states), states_out,
    paste(sprintf("%s=%.0fms", brk$Type, brk$ms), collapse = "  ")
  ))
  invisible(tasks)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: parsec_tasks_to_parquet.r <run_dir> [run_dir ...]")
dbp2xml <- resolve_dbp2xml()
message("dbp2xml: ", dbp2xml)
for (d in args) convert_one(d, dbp2xml)
