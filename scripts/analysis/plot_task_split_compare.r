#!/usr/bin/env Rscript
#
# Comparativo poti vs tupi de ONDE as tarefas executam (contagem CPU vs GPU por
# config), so Cholesky (dpotrf). E o slide que substitui os paineis ST no deck
# final: mostra o under-offload do PaRSEC (device layer estatico) contra o
# dmda/dmdas que enche a GPU -- a perna central da historia "regime decide".
#
# Fonte: application.parquet de cada run. Para o StarPU e o FxT nativo; para o
# PaRSEC a atribuicao de device do application.parquet foi VALIDADA contra os
# contadores nativos do log (bate exato: tupi gd 5913, poti gd 32964/lfq 36539).
# NAO usar states.parquet aqui -- as lanes CUDA dele so tem 18-36% dos kernels
# de GPU (ver memoria parsec-trace-gpu-undercount).
#
#   Rscript scripts/analysis/plot_task_split_compare.r [poti_runs_dir tupi_runs_dir]
#   (default: os dois jobs de trace GPU de 2026-07-01)
#
# Saida: plots/final/task_split_compare_dpotrf.{png,pdf}
#      + plots/final/task_split_compare_dpotrf.csv

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(tidyr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args <- commandArgs(trailingOnly = TRUE)
jobs <- if (length(args) >= 2) {
  c(poti = args[[1]], tupi = args[[2]])
} else {
  c(poti = "data/gpu_traces_poti_20260701_131857/runs",
    tupi = "data/gpu_traces_tupi_20260701_131907/runs")
}
gpu_lbl <- c(poti = "Poti (RTX 4070)", tupi = "Tupi (RTX 4090)")

per_node <- lapply(names(jobs), function(node) {
  node_lbl <- gpu_lbl[[node]]
  runs <- list_runs(jobs[[node]]) %>%
    filter(.data$algo == "dpotrf", .data$application)
  bind_rows(lapply(seq_len(nrow(runs)), function(i) {
    r <- runs[i, ]
    a <- read_parquet(file.path(r$path, "application.parquet"))
    a %>%
      mutate(kernel = norm_kernel(.data$Value)) %>%
      filter(.data$kernel %in% CHOLESKY_KERNELS) %>%
      mutate(class = ifelse(grepl("CUDA", .data$ResourceId), "GPU", "CPU")) %>%
      count(.data$class, name = "tasks") %>%
      mutate(node = node, gpu = node_lbl,
             cfg = paste0(r$runtime, ":", r$scheduler),
             n = r$n, b = r$b)
  }))
})
dat <- bind_rows(per_node) %>%
  group_by(.data$node, .data$cfg) %>%
  mutate(share = .data$tasks / sum(.data$tasks)) %>%
  ungroup() %>%
  mutate(class = factor(.data$class, levels = c("CPU", "GPU")))

print(as.data.frame(dat))
out_dir <- "plots/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(dat, file.path(out_dir, "task_split_compare_dpotrf.csv"),
          row.names = FALSE)

# rotulo do % na GPU acima de cada barra
gpu_pct <- dat %>% filter(.data$class == "GPU") %>%
  group_by(.data$node) %>% mutate(total = .data$tasks / .data$share) %>% ungroup()

# facet com escala livre: o total de tarefas difere ~8x (b=500 vs b=1000)
p <- ggplot(dat, aes(.data$cfg, .data$tasks, fill = .data$class)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = format(.data$tasks, big.mark = ".", decimal.mark = ",")),
            position = position_stack(vjust = 0.5), size = 3.4, colour = "white") +
  geom_text(data = gpu_pct,
            aes(.data$cfg, .data$total,
                label = sprintf("%.0f%% GPU", 100 * .data$share)),
            inherit.aes = FALSE, vjust = -0.4, size = 3.6, fontface = "bold") +
  facet_wrap(~gpu, scales = "free_y") +
  scale_fill_manual(values = c(CPU = "#377eb8", GPU = "#e41a1c")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "tarefas executadas", fill = NULL,
       subtitle = paste("dpotrf FP64, N=40000 (poti b=500, tupi b=1000);",
                        "atribuicao validada contra os contadores nativos")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "task_split_compare_dpotrf.png"), p,
       width = 12, height = 6, dpi = 150)
ggsave(file.path(out_dir, "task_split_compare_dpotrf.pdf"), p,
       width = 12, height = 6, device = grDevices::cairo_pdf)
message("escrito ", out_dir, "/task_split_compare_dpotrf.{png,pdf,csv}")
