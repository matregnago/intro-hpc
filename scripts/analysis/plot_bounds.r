#!/usr/bin/env Rscript
#
# Idea #2 (graficos.md) + extra: how close each execution is to the lower bounds.
# For every run it computes, from the trace:
#
#   makespan  = max(End) - min(Start) over compute states.
#   ABE       = Area Bound Estimate = sum(compute task durations) / #workers.
#               This is exactly what StarVZ's geom_abe computes for a homogeneous
#               (CPU-only) set of resources: the perfect-load-balance lower bound,
#               independent of the DAG. Available for BOTH runtimes from
#               application.parquet alone.
#   CPB       = Critical Path Bound = longest dependency-weighted path through the
#               DAG (infinite-resources lower bound). Needs dag.parquet, so it is
#               drawn only where present (StarPU natively; PaRSEC only where the
#               grapher .dot was collected + parsec_dag_to_parquet.r was run).
#
# Renders plots/bounds_makespan.{png,pdf}: a bar per runtime:scheduler (makespan)
# with the ABE (and CPB) bounds overlaid, annotated with the makespan/ABE ratio
# ("efficiency vs the area bound"). Also prints/writes the numbers.
#
# Usage:  plot_bounds.r [base_dir]
#   defaults to the dpotrf n19200 b480 job (4 runs; StarPU carries dag.parquet).

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(forcats); library(tidyr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (nrow(runs) == 0) stop("no usable runs under ", base_dir)

# makespan + ABE from the unified exec table (ms), both runtimes.
makespan_abe <- function(run_dir) {
  ex <- read_exec(run_dir)
  if (is.null(ex) || nrow(ex) == 0) return(c(makespan = NA, abe = NA, workers = NA))
  workers  <- n_distinct(ex$worker)
  makespan <- max(ex$End) - min(ex$Start)
  abe      <- sum(ex$Duration) / workers
  c(makespan = makespan, abe = abe, workers = workers)
}

# Critical Path Bound: longest dependency-weighted path through dag.parquet.
# StarVZ's dag is bipartite -- timed task nodes plus zero-duration data/other
# nodes (Start/End = NA) that connect them -- so dropping the untimed nodes would
# sever chains. We keep ALL nodes (untimed -> duration 0) and do a Kahn
# topological pass with longest[v] = dur[v] + max over preds u of longest[u].
# Integer-indexed for speed (the n19200 dag is ~49k edges / ~17k nodes).
critical_path <- function(run_dir) {
  f <- file.path(run_dir, "dag.parquet")
  if (!file.exists(f)) return(NA_real_)
  dag <- tryCatch(read_parquet(f), error = function(e) NULL)
  if (is.null(dag) || nrow(dag) == 0) return(NA_real_)

  e <- dag %>%
    transmute(from = as.character(.data$Dependent),
              to   = as.character(.data$JobId)) %>%
    filter(!is.na(.data$from), !is.na(.data$to))
  if (nrow(e) == 0) return(NA_real_)

  # per-task duration (timed nodes); everything else gets 0
  dur_tbl <- dag %>%
    filter(!is.na(.data$Start), !is.na(.data$End)) %>%
    distinct(.data$JobId, .keep_all = TRUE) %>%
    transmute(JobId = as.character(.data$JobId),
              dur = pmax(.data$End - .data$Start, 0))

  nodes <- union(e$from, e$to)
  n     <- length(nodes)
  from_i <- match(e$from, nodes)
  to_i   <- match(e$to,   nodes)
  dur <- numeric(n)
  di <- match(dur_tbl$JobId, nodes)          # timed tasks isolated from any edge -> dropped
  dur[di[!is.na(di)]] <- dur_tbl$dur[!is.na(di)]

  succ  <- split(to_i, from_i)                 # node index -> successor indices
  indeg <- tabulate(to_i, nbins = n)
  longest <- dur
  queue <- integer(n); qh <- 1L; qt <- 0L
  for (v in which(indeg == 0L)) { qt <- qt + 1L; queue[qt] <- v }
  while (qh <= qt) {
    u  <- queue[qh]; qh <- qh + 1L
    su <- succ[[as.character(u)]]
    for (v in su) {
      if (longest[u] + dur[v] > longest[v]) longest[v] <- longest[u] + dur[v]
      indeg[v] <- indeg[v] - 1L
      if (indeg[v] == 0L) { qt <- qt + 1L; queue[qt] <- v }
    }
  }
  if (qt < n) warning("dag not fully ordered (cycle?) in ", run_dir,
                      ": ", qt, "/", n, " nodes")
  max(longest) * time_scale_to_ms(run_dir)     # dag times -> ms (PaRSEC us)
}

res <- bind_rows(lapply(seq_len(nrow(runs)), function(i) {
  r  <- runs[i, ]
  ma <- makespan_abe(r$path)
  tibble(algo = r$algo, cfg = paste0(r$runtime, ":", r$scheduler),
         runtime = r$runtime,
         makespan = ma[["makespan"]], abe = ma[["abe"]],
         workers = ma[["workers"]], cpb = critical_path(r$path),
         gpu = gpu_busy(r$path))
})) %>%
  filter(!is.na(.data$makespan)) %>%
  # ABE (homogeneous area, Sum dur / nworkers) is only a valid bound CPU-only; on
  # GPU runs it can exceed the makespan (work offloaded, but tasks attributed to
  # the CPU threads). There the valid "area" bound is the GPU device-busy floor.
  mutate(is_gpu = !is.na(.data$gpu) & .data$gpu > 0,
         abe = ifelse(.data$is_gpu, NA_real_, .data$abe),
         eff_abe = .data$makespan / .data$abe,
         eff_cpb = .data$makespan / .data$cpb,
         eff_gpu = .data$makespan / .data$gpu)

print(res, n = Inf)
dir.create("plots", showWarnings = FALSE)
write.csv(res, "plots/bounds_summary.csv", row.names = FALSE)

# Bars = makespan; overlay the valid bounds as horizontal markers per bar. ABE is
# NA on GPU and GPU is NA on CPU, so each run shows exactly its valid pair
# (CPU: ABE+CPB; GPU: CPB+GPU) after dropping NA values.
bounds_long <- res %>%
  select(algo = "algo", cfg = "cfg", ABE = "abe", CPB = "cpb", GPU = "gpu") %>%
  pivot_longer(c("ABE", "CPB", "GPU"), names_to = "bound", values_to = "value") %>%
  filter(!is.na(.data$value))

# Annotation: CPB ratio always; GPU-occupancy ratio too where GPU is the area bound.
res <- res %>%
  mutate(lbl = ifelse(.data$is_gpu,
                      sprintf("CPB x%.2f\nGPU x%.2f", .data$eff_cpb, .data$eff_gpu),
                      sprintf("CPB x%.2f", .data$eff_cpb)))

p <- ggplot(res, aes(.data$cfg)) +
  geom_col(aes(y = .data$makespan, fill = .data$runtime), width = 0.6,
           alpha = 0.85) +
  geom_errorbar(data = bounds_long,
                aes(x = .data$cfg, ymin = .data$value, ymax = .data$value,
                    colour = .data$bound), width = 0.6, linewidth = 0.9) +
  geom_text(aes(y = .data$makespan, label = .data$lbl),
            vjust = -0.2, size = 3, lineheight = 0.9) +
  facet_wrap(~algo, scales = "free") +
  scale_colour_manual(values = c(ABE = "black", CPB = "red", GPU = "blue")) +
  labs(title = "Makespan vs limites inferiores (CPB = caminho critico)",
       subtitle = paste0(base_dir,
            "  |  ABE = area (Sum dur/nworkers, so valido na CPU homogenea);  ",
            "GPU = ocupacao do device (piso de throughput, valido na GPU)"),
       x = NULL, y = "tempo (ms)", fill = "runtime", colour = "bound") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p, "bounds_makespan", width = 11, height = 6)
