#!/usr/bin/env Rscript
#
# Build a StarVZ-style tasks.parquet for a PaRSEC run from its grapher .dot
# (structure: dependencies, kernel names, priorities) + binary .prof (timing,
# via dbp2xml), joined on the tpid:tid key both carry. Also writes
# states.parquet: every timed interval in the .prof (compute + runtime/
# scheduler/memory classes) as a long StarVZ-style table.
#
# Iteration is DERIVED (not in the trace): each outer iteration k of a
# right-looking factorization opens with exactly one diagonal kernel, so
# counting those over the DTD submit order recovers k.
#
# Usage:  parsec_tasks_to_parquet.r <run_dir> [run_dir ...]
#   needs cham_*.dot + cham_*.prof-* in run_dir, dbp2xml on PATH or $DBP2XML.

suppressMessages({
  library(xml2); library(arrow); library(dplyr); library(stringr); library(tidyr)
})

# PaRSEC dictionary classes that are not compute kernels; anything not listed
# falls through to "Computation".
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

# Diagonal kernel that opens one outer iteration, by EXACT PaRSEC kernel name:
# Cholesky -> potrf, QR -> geqrt, LU(nopiv) -> getrf_nopiv.
ITER_OPENERS <- c("potrf", "geqrt", "getrf_nopiv")

# Column order/types of StarVZ's StarPU tasks.parquet.
TASKS_COLS <- c(
  "Control", "JobId", "SubmitOrder", "SubmitTime", "MPIRank", "Name", "Model",
  "File", "Line", "Priority", "DependsOn", "DepLabels", "Tag", "WorkerId",
  "MemoryNode", "StartTime", "EndTime", "Footprint", "GFlop", "Iteration",
  "Parameters", "NumaNodes"
)

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

# .prof -> dbp2xml -> one row per timed event across ALL dictionary classes.
parse_prof_events <- function(prof_path, dbp2xml) {
  # dbp2xml hardcodes its output to "out.xml" in the CWD, so run it in a
  # scratch dir with an absolute prof path.
  prof_abs <- normalizePath(prof_path)
  work <- tempfile("dbp2xml_"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old <- setwd(work); on.exit(setwd(old), add = TRUE)
  # dbp2xml carries its own RPATH; the R wrapper's LD_LIBRARY_PATH points at a
  # different glibc that breaks it, so clear it for the child.
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
    # CPU workers identify as "PaRSEC Thread N of VP ...", GPU streams as
    # "GPU <dev>-<stream>". WorkerId stays the CPU thread int (NA for GPU) but
    # ResourceId is always set -- otherwise GPU-executed tasks are dropped.
    cpu <- str_match(ident, "Thread (\\d+)")[, 2]
    gpu <- str_match(ident, "GPU (\\d+)-(\\d+)")
    worker   <- as.integer(cpu)
    resource <- if (!is.na(cpu)) paste0("CPU", cpu) else
                if (!is.na(gpu[1, 1])) paste0("CUDA", gpu[1, 2], "_", gpu[1, 3]) else
                NA_character_
    blocks <- xml_find_all(th, "./KEY")
    if (length(blocks) == 0) return(NULL)
    bind_rows(lapply(blocks, function(b) {
      evs <- xml_find_all(b, "./EVENT")
      if (length(evs) == 0) return(NULL)
      tibble(
        class      = unname(id2name[xml_attr(b, "ID")]),
        key        = xml_text(xml_find_first(evs, "./ID")),
        start_ns   = as.numeric(xml_text(xml_find_first(evs, "./START"))),
        end_ns     = as.numeric(xml_text(xml_find_first(evs, "./END"))),
        WorkerId   = worker,
        ResourceId = resource
      )
    }))
  })
  rows <- bind_rows(rows)
  if (nrow(rows) == 0) {
    return(tibble(class = character(), key = character(), start_ns = numeric(),
                  end_ns = numeric(), WorkerId = integer(), ResourceId = character()))
  }
  rows
}

# Compute kernels only, collapsed to one [start,end] + worker per task key.
# A GPU-executed kernel is profiled TWICE: the CPU manager thread records a
# ~2us submit hook, and the CUDA stream records the real execution span. A
# naive min/max merge would both mis-attribute the task to a CPU lane and
# inflate its span to the whole [CPU submit -> GPU complete] window, so when a
# key has a CUDA-stream event only the CUDA rows are kept.
compute_timing <- function(ev) {
  ev <- ev[event_category(ev$class) == "Computation", , drop = FALSE]
  if (nrow(ev) == 0) {
    return(tibble(key = character(), start_ns = numeric(), end_ns = numeric(),
                  WorkerId = integer(), ResourceId = character()))
  }
  ev <- ev |>
    mutate(is_gpu = !is.na(.data$ResourceId) & startsWith(.data$ResourceId, "CUDA"))
  has_gpu <- ev |> group_by(.data$key) |>
    summarise(gpu = any(.data$is_gpu), .groups = "drop")
  ev |>
    left_join(has_gpu, by = "key") |>
    filter(!.data$gpu | .data$is_gpu) |>
    group_by(.data$key) |>
    summarise(start_ns = min(start_ns), end_ns = max(end_ns),
              WorkerId   = WorkerId[which.min(start_ns)],
              ResourceId = ResourceId[which.min(start_ns)], .groups = "drop")
}

# Every timed interval as a long state row, anchored to t0. Overlapping
# recordings of the same class on a resource are unioned so occupancy is not
# double-counted; genuinely separate events stay distinct rows. Times in us.
build_states <- function(ev, t0) {
  empty <- tibble(Worker = integer(), ResourceId = character(),
                  Type = character(), Value = character(), Start = numeric(),
                  End = numeric(), Duration = numeric(), JobId = character())
  ev <- ev[!is.na(ev$ResourceId), , drop = FALSE]
  if (nrow(ev) == 0) return(empty)
  out <- ev |>
    arrange(ResourceId, class, start_ns, end_ns) |>
    group_by(ResourceId, class) |>
    # new interval whenever this event starts after the running max end so far
    mutate(.grp = cumsum(start_ns > lag(cummax(end_ns), default = -Inf))) |>
    group_by(ResourceId, class, .grp) |>
    summarise(start_ns = min(start_ns), end_ns = max(end_ns),
              WorkerId = first(WorkerId), JobId = first(key), .groups = "drop") |>
    transmute(
      Worker     = as.integer(WorkerId),
      ResourceId = ResourceId,
      Type       = event_category(class),
      Value      = class,
      Start      = (start_ns - t0) / 1000,
      End        = (end_ns   - t0) / 1000,
      Duration   = (end_ns - start_ns) / 1000,
      JobId      = JobId
    ) |>
    arrange(ResourceId, Start) |>
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
    key, WorkerId, ResourceId,
    StartTime = (start_ns - t0) / 1000, EndTime = (end_ns - t0) / 1000
  )

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
    # k = (#openers seen so far) - 1 over submit order; tasks before the first
    # opener (matrix generation) get NA.
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
      ResourceId  = ResourceId,
      MemoryNode  = 0L,
      StartTime   = StartTime,
      EndTime     = EndTime,
      Footprint   = NA_character_,
      GFlop       = NA_real_,
      Iteration   = as.integer(Iteration),
      Parameters  = NA_character_,
      NumaNodes   = NA_character_
    ) |>
    select(all_of(TASKS_COLS), "ResourceId")

  out <- file.path(run_dir, "tasks.parquet")
  write_parquet(tasks, out)

  n_timed <- sum(!is.na(tasks$StartTime))
  niter <- suppressWarnings(max(tasks$Iteration, na.rm = TRUE)) + 1L
  message(sprintf(
    "    tasks=%d (timed=%d, %.0f%%) iters=%s resources=%s prio=[%d..%d] -> %s",
    nrow(tasks), n_timed, 100 * n_timed / nrow(tasks),
    if (is.finite(niter)) niter else "NA",
    paste(sort(unique(na.omit(tasks$ResourceId))), collapse = ","),
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
