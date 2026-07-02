#!/usr/bin/env Rscript
#
# O plot-chave do "por que o parsec:gd ganha no poti/FP64": para cada run de GPU,
# decompoe o trabalho da fatoracao por CLASSE DE RECURSO (CPU vs GPU) e mostra
# quem entregou os flops. Tres visoes + uma tabela:
#
#   gpu_split_contrib.{png,pdf}  -- barras empilhadas: GFlops ENTREGUES por classe
#                                   (gflop da classe / makespan). A soma e o
#                                   desempenho do run; os segmentos mostram o
#                                   trade-off (dmda maximiza GPU, gd maximiza CPU).
#   gpu_split_share.{png,pdf}    -- % dos flops executados na GPU por config, com
#                                   um marcador no split otimo r* daquele run
#                                   (r* = cap_GPU / (cap_GPU + cap_CPU), das taxas
#                                   efetivas medidas no proprio run).
#   gpu_split_rates.{png,pdf}    -- capacidade efetiva por classe: GPU = GFlops
#                                   com o device ativo (uniao de intervalos, imune
#                                   ao overlap dos streams); CPU = n_workers x
#                                   GFlops por worker-segundo ocupado.
#   gpu_split_summary.csv        -- a tabela por run x classe atras das figuras.
#
# Flops por tarefa calculados PELA MESMA formula nos dois runtimes (kernel x b do
# nome do run: gemm 2b^3, trsm/syrk b^3, potrf b^3/3, getrf_nopiv 2b^3/3), em vez
# da coluna GFlop do StarVZ (que o PaRSEC nao tem) -- tratamento identico, sem
# vies. Kernels de init (plgsy etc.) ficam fora de tudo (rodam antes da regiao
# cronometrada).
#
# Uso:  plot_gpu_split.r [base_dir]
# base_dir default: o job de traces GPU do poti (N=40000 b=500, 8 runs).

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(tidyr); library(stringr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/gpu_traces_poti_20260701_131857/runs"

# coeficiente de flops em unidades de b^3 (FP64/FP32 indiferente: e contagem)
FLOP_COEF <- c(gemm = 2, trsm = 1, syrk = 1, herk = 1, trmm = 1, potrf = 1 / 3,
               getrf_nopiv = 2 / 3, getrf = 2 / 3)

runs <- list_runs(base_dir) %>% filter(.data$application | .data$tasks)
if (nrow(runs) == 0) stop("no runs with application/tasks parquet under ", base_dir)
message("runs: ", paste(runs$dir, collapse = ", "))

per_run <- lapply(seq_len(nrow(runs)), function(i) {
  r  <- runs[i, ]
  ex <- read_exec(r$path)
  if (is.null(ex)) return(NULL)
  ex <- ex %>%
    filter(!.data$kernel %in% INIT_KERNELS,
           .data$kernel %in% names(FLOP_COEF)) %>%
    mutate(class = ifelse(grepl("CUDA", .data$worker), "GPU", "CPU"),
           gflop = FLOP_COEF[.data$kernel] * as.numeric(r$b)^3 / 1e9)
  if (nrow(ex) == 0) return(NULL)
  makespan_s <- (max(ex$End) - min(ex$Start)) / 1000
  ex %>%
    group_by(.data$class) %>%
    summarise(tasks    = n(),
              workers  = n_distinct(.data$worker),
              busy_s   = sum(.data$Duration) / 1000,
              active_s = union_length(.data$Start, .data$End) / 1000,
              gflop    = sum(.data$gflop), .groups = "drop") %>%
    mutate(algo = r$algo, runtime = r$runtime,
           cfg = paste0(r$runtime, ":", r$scheduler),
           makespan_s = makespan_s)
})
dat <- bind_rows(per_run)
if (nrow(dat) == 0) stop("no compute spans found under ", base_dir)

dat <- dat %>%
  mutate(
    # GFlops entregues por classe ao longo do run inteiro
    contrib_gfps = .data$gflop / .data$makespan_s,
    # capacidade efetiva: GPU pela uniao (streams overlapam); CPU por
    # worker-segundo ocupado, escalado pelo n de workers da classe
    capacity_gfps = ifelse(.data$class == "GPU",
                           .data$gflop / .data$active_s,
                           .data$workers * .data$gflop / .data$busy_s)
  ) %>%
  group_by(.data$algo, .data$cfg) %>%
  mutate(share_flops = .data$gflop / sum(.data$gflop)) %>%
  ungroup()

# split otimo por run: fracao da capacidade total que esta na GPU
rstar <- dat %>%
  select("algo", "cfg", "class", "capacity_gfps") %>%
  pivot_wider(names_from = "class", values_from = "capacity_gfps") %>%
  mutate(rstar = .data$GPU / (.data$GPU + .data$CPU)) %>%
  select("algo", "cfg", "rstar")

dat <- dat %>% left_join(rstar, by = c("algo", "cfg"))
print(as.data.frame(dat), digits = 4)
write.csv(dat, file.path(plots_dir(), "gpu_split_summary.csv"), row.names = FALSE)

cls_fill <- scale_fill_manual(values = c(CPU = "#377eb8", GPU = "#e41a1c"))

# ---- (1) GFlops entregues, empilhado por classe ----
tot <- dat %>% group_by(.data$algo, .data$cfg) %>%
  summarise(total = sum(.data$contrib_gfps), .groups = "drop")
p1 <- ggplot(dat, aes(.data$cfg, .data$contrib_gfps, fill = .data$class)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f", .data$contrib_gfps)),
            position = position_stack(vjust = 0.5), size = 3.2, colour = "white") +
  geom_text(data = tot, aes(.data$cfg, .data$total,
                            label = sprintf("%.0f", .data$total)),
            inherit.aes = FALSE, vjust = -0.4, size = 3.4, fontface = "bold") +
  facet_wrap(~algo) +
  cls_fill +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(title = "Quem entrega os GFlops: contribuicao de CPU e GPU por config",
       subtitle = paste0(base_dir, "  |  total = desempenho do run; ",
                         "dmda ganha na GPU mas perde mais na CPU"),
       x = NULL, y = "GFlops entregues (gflops da classe / makespan)",
       fill = "classe") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
save_plot(p1, "gpu_split_contrib", width = 11, height = 6)

# ---- (2) % dos flops na GPU vs split otimo r* ----
shr <- dat %>% filter(.data$class == "GPU")
p2 <- ggplot(shr, aes(.data$cfg, 100 * .data$share_flops)) +
  geom_col(width = 0.6, alpha = 0.9, fill = "#e41a1c") +
  geom_point(aes(y = 100 * .data$rstar), shape = 18, size = 4, colour = "black") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * .data$share_flops)),
            vjust = -0.5, size = 3.3) +
  facet_wrap(~algo) +
  labs(title = "Fracao dos flops executada na GPU vs split otimo",
       subtitle = paste0(base_dir,
                         "  |  losango = r* = cap_GPU/(cap_GPU+cap_CPU) medido no run"),
       caption = paste("acima do losango = GPU sobrecarregada (e o gargalo,",
                       "eff_gpu ~1.0 em plot_bounds) -- quem carrega menos a GPU",
                       "termina antes"),
       x = NULL, y = "% dos GFlops na GPU") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p2, "gpu_split_share", width = 11, height = 6)

# ---- (3) capacidade efetiva por classe ----
p3 <- ggplot(dat, aes(.data$cfg, .data$capacity_gfps, fill = .data$class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f", .data$capacity_gfps)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3.1) +
  facet_wrap(~algo) +
  cls_fill +
  labs(title = "Capacidade efetiva por classe de recurso",
       subtitle = paste0(base_dir, "  |  GPU = GFlops c/ device ativo; ",
                         "CPU = workers x taxa por worker ocupado"),
       caption = paste("a CPU do StarPU rende menos que a do PaRSEC (4 threads",
                       "dedicadas aos workers CUDA oversubscrevem os cores);",
                       "a GPU do StarPU rende mais (prefetch + so 3.4% de re-stage)"),
       x = NULL, y = "GFlops", fill = "classe") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
save_plot(p3, "gpu_split_rates", width = 11, height = 6)
