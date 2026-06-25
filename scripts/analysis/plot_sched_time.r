#!/usr/bin/env Rscript
#
# Per-task scheduling delay: the gap between when a task became *ready* (its last
# dependency finished) and when it actually *started*:
#
#   sched(T) = StartTime(T) - max_{D in deps(T)} EndTime(D)
#
# This isolates the runtime's scheduling/queueing latency from compute time. It
# reads tasks.parquet (DependsOn + StartTime/EndTime), which exists for StarPU
# (StarVZ) and for PaRSEC (scripts/analysis/parsec_tasks_to_parquet.r), so the
# same metric is computed identically for both runtimes.
#
# Run dirs are hardcoded below (edit when pointing at a new job), mirroring
# plot_traces.r. Outputs plots/sched_time.{png,pdf}.

suppressMessages({
  library(arrow); library(dplyr); library(tidyr); library(stringr); library(ggplot2)
})

runs <- c(
  "parsec / spotrf (n3840)"  = "data/parsec_dag/runs/spotrf_lfq_n3840_b960",
  "parsec / sgeqrf (n3840)"  = "data/parsec_dag/runs/sgeqrf_lfq_n3840_b960",
  "starpu / dpotrf (n19200)" = "data/parcial/chameleon-runtimes_782088/runs/0001_starpu_lws_dpotrf_n19200_b480_rep1"
)

sched_of <- function(run_dir) {
  t <- read_parquet(file.path(run_dir, "tasks.parquet"))
  end_of <- t |> select(JobId, EndTime)              # dependency -> its end time
  t |>
    filter(!is.na(StartTime), !is.na(DependsOn)) |>
    select(JobId, Name, StartTime, DependsOn) |>
    separate_rows(DependsOn, sep = " ") |>
    left_join(end_of, by = c("DependsOn" = "JobId")) |>
    group_by(JobId, Name, StartTime) |>
    summarise(dep_ready = suppressWarnings(max(EndTime, na.rm = TRUE)), .groups = "drop") |>
    filter(is.finite(dep_ready)) |>                  # keep tasks with >=1 timed dep
    mutate(sched = StartTime - dep_ready)
}

all <- bind_rows(lapply(names(runs), function(lbl) {
  sched_of(runs[[lbl]]) |> mutate(Run = lbl)
})) |>
  mutate(Run = factor(Run, levels = names(runs)))

# numeric summary to stdout
cat("\n=== scheduling delay (us) por run/kernel ===\n")
print(all |>
  group_by(Run, Name) |>
  summarise(n = n(), median = round(median(sched), 1),
            mean = round(mean(sched), 1), max = round(max(sched), 1),
            negativos = sum(sched < 0), .groups = "drop") |>
  as.data.frame())

# per-task scheduling delay, grouped by kernel, one facet per run
p <- ggplot(all, aes(x = Name, y = sched, color = Name)) +
  geom_jitter(width = 0.25, alpha = 0.35, size = 0.8) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  facet_wrap(~Run, scales = "free", nrow = 1) +
  labs(x = "kernel", y = "tempo de escalonamento  StartTime - max(End dep)  [us]",
       title = "Tempo de escalonamento por tarefa (atraso entre ficar pronta e iniciar)") +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 30, hjust = 1))

dir.create("plots", showWarnings = FALSE)
ggsave("plots/sched_time.png", p, width = 14, height = 5, dpi = 300)
ggsave("plots/sched_time.pdf", p, width = 14, height = 5, device = cairo_pdf)
cat("\nsalvo: plots/sched_time.png e .pdf\n")
