#!/usr/bin/env Rscript
#
# Plot the reconstructed "Ready" scheduler counter for a PaRSEC run, two ways:
#   1. plots/ready_panel_<slug>.{png,pdf} -- the native StarVZ panel_ready, fed
#      with the variable.parquet reconstructed by parsec_ready_to_parquet.r.
#   2. plots/ready_raw_<slug>.{png,pdf}   -- a sanity plot of the raw step
#      function: Ready (eligible-but-not-running) vs Running (executing), with the
#      worker count as a reference line.
#
# Ready is reconstructed on the fly from the run's tasks.parquet (DAG + timing);
# nothing is written into data/. See parsec_ready_to_parquet.r for the method.
#
# Usage:
#   plot_ready.r [run_dir] [label]
# Defaults to the manual_traces lfq spotrf n19200 run.

suppressMessages({
  library(arrow); library(dplyr); library(tibble)
  library(ggplot2); library(starvz)
})

# Pull in reconstruct_ready() without triggering its CLI.
repo <- tryCatch(dirname(dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))))),
  error = function(e) ".")
source(file.path(repo, "scripts/analysis/parsec_ready_to_parquet.r"))

args     <- commandArgs(trailingOnly = TRUE)
run_dir  <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs/0005_parsec_lfq_spotrf_n19200_b960_rep1"
label    <- if (length(args) >= 2) args[[2]] else "PaRSEC lfq spotrf n19200 b960"
slug     <- gsub("[^A-Za-z0-9]+", "_", label)

dir.create("plots", showWarnings = FALSE)

tk <- read_parquet(file.path(run_dir, "tasks.parquet")) %>%
  filter(!is.na(.data$StartTime), !is.na(.data$EndTime))
nworkers <- n_distinct(tk$WorkerId)
gmin     <- min(tk$StartTime)
makespan <- max(tk$EndTime) - gmin

Variable <- reconstruct_ready(tk)

# ---- (1) native StarVZ panel_ready ----
svz <- list(
  Variable = Variable,
  config = list(base_size = 14, expand = 0.05,
                ready = list(step = makespan / 400, legend = TRUE,
                             limit = NA, lack_ready = list(active = FALSE)))
)
p_panel <- panel_ready(svz) + ggtitle(paste0("panel_ready (StarVZ) -- ", label))
ggsave(paste0("plots/ready_panel_", slug, ".png"), p_panel,
       width = 12, height = 4.5, dpi = 150, limitsize = FALSE)
ggsave(paste0("plots/ready_panel_", slug, ".pdf"), p_panel,
       width = 12, height = 4.5, device = cairo_pdf)

# ---- (2) raw sanity plot: Ready vs Running ----
running <- bind_rows(
  tibble(t = tk$StartTime - gmin, d =  1L),
  tibble(t = tk$EndTime   - gmin, d = -1L)
) %>%
  group_by(.data$t) %>% summarise(d = sum(.data$d), .groups = "drop") %>%
  arrange(.data$t) %>% mutate(run = cumsum(.data$d))

p_raw <- ggplot() +
  geom_step(data = Variable, aes(.data$Start, .data$Value,
                                 colour = "Ready (reconstruido)")) +
  geom_step(data = running, aes(.data$t, .data$run,
                                colour = "Running (em execucao)")) +
  geom_hline(yintercept = nworkers, linetype = "dashed") +
  annotate("text", x = 0, y = nworkers, vjust = -0.4, hjust = 0,
           label = paste0(nworkers, " workers"), size = 3) +
  labs(title = paste0(label, " -- Ready reconstruido do DAG + timing"),
       x = "tempo (unidades do trace, t0 = 1o start)",
       y = "no de tarefas", colour = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "top")
ggsave(paste0("plots/ready_raw_", slug, ".png"), p_raw,
       width = 12, height = 5, dpi = 150, limitsize = FALSE)
ggsave(paste0("plots/ready_raw_", slug, ".pdf"), p_raw,
       width = 12, height = 5, device = cairo_pdf)

message(sprintf("salvo: plots/ready_{panel,raw}_%s.{png,pdf}  (Ready max=%g mean=%.1f, %d workers)",
                slug, max(Variable$Value),
                weighted.mean(Variable$Value, Variable$Duration), nworkers))
