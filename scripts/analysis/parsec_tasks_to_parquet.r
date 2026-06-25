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
#   - Model/File/Line/Footprint/GFlop/Parameters/SubmitTime/Tag/Iteration/
#     NumaNodes/Control               <- not in the PaRSEC trace -> NA
#       (GFlop is derivable analytically per kernel+tile; left NA for now.)
#
# Usage:  parsec_tasks_to_parquet.r <run_dir> [run_dir ...]
#   needs cham_*.dot + cham_*.prof-* in run_dir, dbp2xml on PATH or $DBP2XML.
#   Writes <run_dir>/tasks.parquet. See parsec_dag_to_parquet.r for the
#   companion dag.parquet and the shared dbp2xml/out.xml gotcha.

suppressMessages({
  library(xml2); library(arrow); library(dplyr); library(stringr); library(tidyr)
})

INTERNAL_CLASSES <- c(
  "MEMALLOC", "Sched POLL", "Sched PUSH", "Sched SLEEP", "Queue REMOVE",
  "ARENA_MEMORY", "ARENA_ACTIVE_SET", "TASK_MEMORY", "Device delegate",
  "RELEASE_DEPS", "ACTIVATE_CB", "DATA_FLUSH"
)

# Final schema (column order + types) of StarVZ's StarPU tasks.parquet.
TASKS_COLS <- c(
  "Control", "JobId", "SubmitOrder", "SubmitTime", "MPIRank", "Name", "Model",
  "File", "Line", "Priority", "DependsOn", "DepLabels", "Tag", "WorkerId",
  "MemoryNode", "StartTime", "EndTime", "Footprint", "GFlop", "Iteration",
  "Parameters", "NumaNodes"
)

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

# .prof -> dbp2xml -> per-task span + worker (tpid:tid -> start/end ns, WorkerId)
parse_prof <- function(prof_path, dbp2xml) {
  prof_abs <- normalizePath(prof_path)
  work <- tempfile("dbp2xml_"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old <- setwd(work); on.exit(setwd(old), add = TRUE)
  system2(dbp2xml, shQuote(prof_abs), stdout = FALSE, stderr = FALSE)
  xml_path <- file.path(work, "out.xml")
  if (!file.exists(xml_path) || file.info(xml_path)$size == 0) {
    stop("dbp2xml produced no out.xml for ", prof_path)
  }
  doc <- read_xml(xml_path)

  keys <- xml_find_all(doc, "//DICTIONARY/KEY")
  dict <- tibble(id = xml_attr(keys, "ID"),
                 name = xml_text(xml_find_first(keys, "./NAME")))
  compute_ids <- dict$id[!(dict$name %in% INTERNAL_CLASSES)]

  rows <- lapply(xml_find_all(doc, "//THREAD"), function(th) {
    ident  <- xml_text(xml_find_first(th, "./IDENTIFIER"))
    worker <- as.integer(str_match(ident, "Thread (\\d+)")[, 2])
    blocks <- xml_find_all(th, "./KEY")
    blocks <- blocks[xml_attr(blocks, "ID") %in% compute_ids]
    if (length(blocks) == 0) return(NULL)
    bind_rows(lapply(blocks, function(b) {
      evs <- xml_find_all(b, "./EVENT")
      if (length(evs) == 0) return(NULL)
      tibble(
        key      = xml_text(xml_find_first(evs, "./ID")),
        start_ns = as.numeric(xml_text(xml_find_first(evs, "./START"))),
        end_ns   = as.numeric(xml_text(xml_find_first(evs, "./END"))),
        WorkerId = worker
      )
    }))
  })
  rows <- bind_rows(rows)
  if (nrow(rows) == 0) {
    return(tibble(key = character(), start_ns = numeric(),
                  end_ns = numeric(), WorkerId = integer()))
  }
  rows |>
    group_by(key) |>
    summarise(start_ns = min(start_ns), end_ns = max(end_ns),
              WorkerId = WorkerId[which.min(start_ns)], .groups = "drop")
}

convert_one <- function(run_dir, dbp2xml) {
  dot  <- list.files(run_dir, pattern = "\\.dot$", full.names = TRUE)
  prof <- list.files(run_dir, pattern = "\\.prof-", full.names = TRUE)
  prof <- prof[!grepl("\\.tar\\.gz$", prof)]
  if (length(dot) == 0)  { message("  no .dot in ", run_dir, ", skip");  return(invisible()) }
  if (length(prof) == 0) { message("  no .prof in ", run_dir, ", skip"); return(invisible()) }

  message("==> ", run_dir)
  g <- parse_dot(dot[1])
  t <- parse_prof(prof[1], dbp2xml)

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
      Iteration   = NA_integer_,
      Parameters  = NA_character_,
      NumaNodes   = NA_character_
    ) |>
    select(all_of(TASKS_COLS))

  out <- file.path(run_dir, "tasks.parquet")
  write_parquet(tasks, out)

  n_timed <- sum(!is.na(tasks$StartTime))
  message(sprintf(
    "    tasks=%d (timed=%d, %.0f%%) workers=%s prio=[%d..%d] -> %s",
    nrow(tasks), n_timed, 100 * n_timed / nrow(tasks),
    paste(sort(unique(na.omit(tasks$WorkerId))), collapse = ","),
    min(tasks$Priority, na.rm = TRUE), max(tasks$Priority, na.rm = TRUE), out
  ))
  invisible(tasks)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: parsec_tasks_to_parquet.r <run_dir> [run_dir ...]")
dbp2xml <- resolve_dbp2xml()
message("dbp2xml: ", dbp2xml)
for (d in args) convert_one(d, dbp2xml)
