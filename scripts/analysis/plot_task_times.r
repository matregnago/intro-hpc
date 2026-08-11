#!/usr/bin/env Rscript
#
# Duracao media por tipo de tarefa, comparada entre runtimes/escalonadores, a
# partir do application.parquet de cada run (fallback: tasks.parquet -- ver
# read_exec). Usado nos slides, nao no artigo.
#
# Uso:  plot_task_times.r [base_dir] [algo]
# Saida: <PLOTS_DIR>/task_times_mean.{png,pdf} + task_times_summary.csv

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
source(file.path(script_dir, "plot_style.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"
algo_sel <- if (length(args) >= 2) args[[2]] else NA

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (!is.na(algo_sel)) runs <- runs %>% filter(.data$algo == algo_sel)
if (nrow(runs) == 0) stop("no runs with application/tasks parquet under ", base_dir,
                          if (!is.na(algo_sel)) paste0(" (algo=", algo_sel, ")") else "")
message("runs: ", paste(runs$dir, collapse = ", "))

read_kernels <- function(row) {
  ex <- read_exec(row$path)
  if (is.null(ex)) return(NULL)
  ex %>%
    transmute(.data$kernel, .data$Duration, .data$src,
              class = ifelse(grepl("CUDA", .data$worker), "GPU", "CPU"),
              runtime = row$runtime, scheduler = row$scheduler,
              algo = row$algo, cfg = cfg_label(row$runtime, row$scheduler))
}

dat <- bind_rows(lapply(seq_len(nrow(runs)), function(i) read_kernels(runs[i, ])))
if (nrow(dat) == 0) stop("no compute-kernel states found")
message("fontes: ", paste(unique(paste0(dat$runtime, "=", dat$src)), collapse = ", "))

core <- dat %>% filter(!.data$kernel %in% INIT_KERNELS)

# Os spans de GPU do PaRSEC sao abertos/fechados no host e se sobrepoem entre as
# streams, inflando a duracao (~2.9x observado). Compara-los com a GPU do StarPU
# sugeriria que o mesmo cuBLAS dgemm e 2.9x mais lento num runtime, o que e
# falso. Ficam no CSV como registro bruto, mas fora da figura.
n_drop <- sum(core$runtime == "parsec" & core$class == "GPU")
if (n_drop > 0)
  message("suprimindo ", n_drop, " spans de GPU do PaRSEC da figura; seguem no CSV")
core <- core %>% filter(!(.data$runtime == "parsec" & .data$class == "GPU"))

summary_tbl <- dat %>%
  group_by(.data$algo, .data$cfg, .data$class, .data$kernel) %>%
  summarise(n = n(), mean_ms = mean(.data$Duration),
            median_ms = median(.data$Duration), sd_ms = sd(.data$Duration),
            .groups = "drop") %>%
  arrange(.data$algo, .data$kernel, .data$cfg, .data$class)
print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(plots_dir(), "task_times_summary.csv"), row.names = FALSE)

# faceta = algo + kernel; com um unico algoritmo, so o kernel
one_algo <- n_distinct(core$algo) == 1
core <- core %>%
  mutate(panel = if (one_algo) .data$kernel
         else paste0(.data$algo, ": ", .data$kernel))

message("para o caption: ", base_dir, " | classe = onde o span executou; ",
        "cross-runtime so CPU vs CPU")

mean_tbl <- core %>%
  group_by(.data$panel, .data$cfg, .data$class) %>%
  summarise(mean_us = mean(.data$Duration),
            se = sd(.data$Duration) / sqrt(n()), .groups = "drop")

p_mean <- ggplot(mean_tbl, aes(.data$cfg, .data$mean_us, fill = .data$class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(aes(ymin = .data$mean_us - .data$se,
                    ymax = .data$mean_us + .data$se),
                position = position_dodge(width = 0.8), width = 0.3) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "duracao media (ms)") +
  theme_sscad(legend = "bottom") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_mean, "task_times_mean", width = 12, height = 8)
