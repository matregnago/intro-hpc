#!/usr/bin/env Rscript
#
# Convert a PaRSEC run's task-graph .dot + binary .prof into a StarVZ-style
# dag.parquet (columns: JobId, Dependent, Start, End, Cost, Value), the same
# table StarVZ builds from StarPU's dag.dot. With it, the DAG-based analyses
# StarVZ does for StarPU (critical path, gantt overlay, ...) also run on PaRSEC.
#
# Why two inputs:
#   - the .dot (PaRSEC grapher, built with PARSEC_PROF_GRAPHER=ON) carries the
#     graph STRUCTURE + kernel name. Each node's tooltip is
#     "tpid=<p>:did=<d>:tname=<kernel>:tid=<t>" and edges are dependencies.
#   - the .prof carries the per-task TIMING. `dbp2xml` dumps it to XML where each
#     <EVENT> under a compute-class <KEY> has <ID>tpid:tid</ID> + <START>/<END>
#     in ns. tpid:tid is the SAME key the .dot uses -> we join timing onto
#     structure. (dbp2paje drops this id, which is why the starvz application
#     trace has JobId = NA and can't be joined directly.)
#
# StarVZ conventions reproduced here:
#   - one row per dependency edge: JobId = task, Dependent = its predecessor.
#   - Start/End/Value describe the JobId task; times in microseconds, shifted so
#     the earliest profiled event is 0.
#   - Cost = Start - End (negative duration; used as edge weight so a
#     shortest-path on Cost is the timing critical path). 0 when untimed.
#   - tasks with no profiling event (init splgsy/splrnt, parsec_dtd_data_flush)
#     keep Start/End = NA, Cost = 0, mirroring StarPU's data nodes.
#
# Usage:
#   parsec_dag_to_parquet.r <run_dir> [run_dir ...]
# Each run_dir must contain cham_*.dot and cham_*.prof-*; dbp2xml must be on
# PATH or pointed to by $DBP2XML. Writes <run_dir>/dag.parquet.

suppressMessages({
  library(xml2)
  library(arrow)
  library(dplyr)
  library(stringr)
})

# PaRSEC profiling dictionary classes that are runtime bookkeeping, not tasks.
# Everything else in the dictionary is a task class (potrf, gemm, geqrt, ...).
INTERNAL_CLASSES <- c(
  "MEMALLOC", "Sched POLL", "Sched PUSH", "Sched SLEEP", "Queue REMOVE",
  "ARENA_MEMORY", "ARENA_ACTIVE_SET", "TASK_MEMORY", "Device delegate",
  "RELEASE_DEPS", "ACTIVATE_CB", "DATA_FLUSH"
)

resolve_dbp2xml <- function() {
  env <- Sys.getenv("DBP2XML")
  if (nzchar(env) && file.exists(env)) return(env)
  found <- Sys.which("dbp2xml")
  if (nzchar(found)) return(unname(found))
  # Fallback: build the flake's parsec package and use its dbp2xml (mirrors
  # scripts/analysis/trace_phase1.sh). --impure: parsec.nix reads $PWD.
  message("dbp2xml not on PATH; building .#parsec (first run only)...")
  out <- suppressWarnings(system2(
    "nix",
    c("build", "--impure", "--no-link", "--print-out-paths", ".#parsec"),
    stdout = TRUE, stderr = FALSE
  ))
  out <- tail(out[nzchar(out)], 1)
  bin <- if (length(out)) file.path(out, "bin", "dbp2xml") else ""
  if (nzchar(bin) && file.exists(bin)) return(bin)
  stop("dbp2xml not found: put it on PATH, set $DBP2XML, or run from a tree ",
       "where `nix build --impure .#parsec` works")
}

# --- .dot: nodes (tpid:tid -> kernel) and edges (predecessor -> task) ----------
parse_dot <- function(dot_path) {
  ln <- readLines(dot_path, warn = FALSE)

  # node declarations: <name> [ ... tooltip="tpid=P:did=D:tname=K:tid=T" ... ];
  decl <- str_match(
    ln, '^\\s*(\\S+)\\s+\\[.*tooltip="tpid=(\\d+):did=\\d+:tname=([^:"]+):tid=(\\d+)"'
  )
  decl <- decl[!is.na(decl[, 1]), , drop = FALSE]
  nodes <- tibble(
    name  = decl[, 2],
    key   = paste(decl[, 3], decl[, 5], sep = ":"),  # tpid:tid
    tname = decl[, 4]
  ) |> distinct(name, .keep_all = TRUE)
  name2key   <- setNames(nodes$key, nodes$name)
  key2tname  <- setNames(nodes$tname, nodes$key)

  # edges: <src> -> <dst>
  em <- str_match(ln, '^\\s*(\\S+)\\s*->\\s*(\\S+)')
  em <- em[!is.na(em[, 1]), , drop = FALSE]
  edges <- tibble(
    src_key = name2key[em[, 2]],
    dst_key = name2key[em[, 3]]
  ) |> filter(!is.na(src_key), !is.na(dst_key))

  list(nodes = nodes, edges = edges, key2tname = key2tname)
}

# --- .prof -> dbp2xml -> per-task timing (tpid:tid -> Start/End in ns) ----------
parse_prof_timing <- function(prof_path, dbp2xml) {
  # dbp2xml ignores stdout and hardcodes its output to "out.xml" in the CWD
  # (tools/profiling/dbp2xml.c: dump_xml("out.xml", dbp)). Run it in a scratch
  # dir, with an absolute prof path, and read that out.xml back.
  prof_abs <- normalizePath(prof_path)
  work <- tempfile("dbp2xml_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE)

  system2(dbp2xml, shQuote(prof_abs), stdout = FALSE, stderr = FALSE)
  xml_path <- file.path(work, "out.xml")
  if (!file.exists(xml_path) || file.info(xml_path)$size == 0) {
    stop("dbp2xml produced no out.xml for ", prof_path)
  }
  doc <- read_xml(xml_path)

  # dictionary: KEY ID -> class name; keep only task (non-internal) classes
  keys <- xml_find_all(doc, "//DICTIONARY/KEY")
  dict <- tibble(
    id   = xml_attr(keys, "ID"),
    name = xml_text(xml_find_first(keys, "./NAME"))
  )
  compute_ids <- dict$id[!(dict$name %in% INTERNAL_CLASSES)]

  # each <THREAD>/<KEY ID=k> block holds the <EVENT>s of class k for one worker.
  # Pull (key=tpid:tid, START, END) per block, tagging it with its class id.
  blocks <- xml_find_all(doc, "//THREAD/KEY")
  ev <- lapply(blocks, function(b) {
    evs <- xml_find_all(b, "./EVENT")
    if (length(evs) == 0) return(NULL)
    tibble(
      class_id = xml_attr(b, "ID"),
      key      = xml_text(xml_find_first(evs, "./ID")),   # tpid:tid
      start_ns = as.numeric(xml_text(xml_find_first(evs, "./START"))),
      end_ns   = as.numeric(xml_text(xml_find_first(evs, "./END")))
    )
  })
  ev <- bind_rows(ev)
  if (nrow(ev) == 0) return(tibble(key = character(), start_ns = numeric(), end_ns = numeric()))
  ev |>
    filter(class_id %in% compute_ids) |>
    group_by(key) |>                                      # one span per task
    summarise(start_ns = min(start_ns), end_ns = max(end_ns), .groups = "drop")
}

convert_one <- function(run_dir, dbp2xml) {
  dot  <- list.files(run_dir, pattern = "\\.dot$", full.names = TRUE)
  prof <- list.files(run_dir, pattern = "\\.prof-", full.names = TRUE)
  prof <- prof[!grepl("\\.tar\\.gz$", prof)]
  if (length(dot) == 0)  { message("  no .dot in ", run_dir, ", skip");  return(invisible()) }
  if (length(prof) == 0) { message("  no .prof in ", run_dir, ", skip"); return(invisible()) }

  message("==> ", run_dir)
  g <- parse_dot(dot[1])
  t <- parse_prof_timing(prof[1], dbp2xml)

  # ns -> us, shifted so the earliest profiled event is 0 (StarVZ origin)
  t0 <- if (nrow(t)) min(t$start_ns) else 0
  timing <- t |>
    transmute(key, Start = (start_ns - t0) / 1000, End = (end_ns - t0) / 1000)
  start_of <- setNames(timing$Start, timing$key)
  end_of   <- setNames(timing$End,   timing$key)

  # column order/types mirror StarVZ's StarPU dag.parquet so starvz_read() and
  # the DAG panels consume it unchanged.
  dag <- g$edges |>
    transmute(
      JobId     = as.character(dst_key),
      Dependent = as.factor(src_key),
      Start     = unname(start_of[dst_key]),
      End       = unname(end_of[dst_key]),
      Cost      = ifelse(is.na(unname(start_of[dst_key])), 0,
                         unname(start_of[dst_key]) - unname(end_of[dst_key])),
      Value     = as.factor(unname(g$key2tname[dst_key]))
    )

  out <- file.path(run_dir, "dag.parquet")
  write_parquet(dag, out)

  n_tasks <- length(unique(g$nodes$key))
  n_timed <- nrow(timing)
  n_job   <- length(unique(dag$JobId))
  hit     <- mean(!is.na(dag$Start))
  message(sprintf(
    "    nodes=%d edges=%d timed_tasks=%d | dag rows=%d distinct JobId=%d timed=%.0f%% -> %s",
    n_tasks, nrow(g$edges), n_timed, nrow(dag), n_job, 100 * hit, out
  ))
  invisible(dag)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("usage: parsec_dag_to_parquet.r <run_dir> [run_dir ...]")
}
dbp2xml <- resolve_dbp2xml()
message("dbp2xml: ", dbp2xml)
for (d in args) convert_one(d, dbp2xml)
