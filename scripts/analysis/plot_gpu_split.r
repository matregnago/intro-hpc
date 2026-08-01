#!/usr/bin/env Rscript
#
# O plot-chave do "por que o parsec:gd ganha no poti/FP64": para cada run de GPU,
# decompoe o trabalho da fatoracao por CLASSE DE RECURSO (CPU vs GPU) e mostra
# quem entregou os flops. Tres visoes + uma tabela:
#
#   gpu_split_contrib.{png,pdf}  -- barras empilhadas: GFLOPS ENTREGUES por classe
#                                   (gflop da classe / makespan). A soma e o
#                                   desempenho do run; os segmentos mostram o
#                                   trade-off (dmda maximiza GPU, gd maximiza CPU).
#                                   E a Figura 3 do artigo.
#   gpu_split_share.{png,pdf}    -- % dos flops executados na GPU por config, com
#                                   um marcador no split otimo r* daquele run
#                                   (r* = cap_GPU / (cap_GPU + cap_CPU), das taxas
#                                   efetivas medidas no proprio run).
#   gpu_split_rates.{png,pdf}    -- capacidade efetiva por classe: GPU = GFLOPS
#                                   com o device ativo (uniao de intervalos, imune
#                                   ao overlap dos streams); CPU = n_workers x
#                                   GFLOPS por worker-segundo ocupado.
#   gpu_split_summary.csv        -- a tabela por run x classe atras das figuras.
#
# Flops por tarefa calculados PELA MESMA formula nos dois runtimes (kernel x b do
# nome do run: gemm 2b^3, trsm/syrk b^3, potrf b^3/3, getrf_nopiv 2b^3/3), em vez
# da coluna GFlop do StarVZ (que o PaRSEC nao tem) -- tratamento identico, sem
# vies. Kernels de init (plgsy etc.) ficam fora de tudo (rodam antes da regiao
# cronometrada).
#
# Sem titulo/subtitulo dentro dos graficos (pedido do revisor): o que era titulo
# vai para o stderr, para ser colado no caption do artigo. Estilo (paleta, fonte,
# ordem das configs, nomes dos algoritmos) vem todo de plot_style.r.
#
# Uso:  plot_gpu_split.r [base_dir]
#   base_dir com run dirs  -> recalcula tudo dos rastros e reescreve o summary.
#   base_dir com gpu_split_summary.csv (e sem run dirs) -> apenas replota daquele
#   CSV, o que permite regerar a Figura 3 sem os dados brutos do PCAD.
# base_dir default: o job de traces GPU do poti (N=40000 b=500, 8 runs).

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2); library(tidyr); library(stringr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
source(file.path(script_dir, "plot_style.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/gpu_traces_poti_20260701_131857/runs"

# coeficiente de flops em unidades de b^3 (FP64/FP32 indiferente: e contagem)
FLOP_COEF <- c(gemm = 2, trsm = 1, syrk = 1, herk = 1, trmm = 1, potrf = 1 / 3,
               getrf_nopiv = 2 / 3, getrf = 2 / 3)

#' Recalcula a decomposicao a partir dos rastros de cada run.
from_traces <- function(base_dir) {
  # list_runs() aborta quando o diretorio nao tem run dirs; aqui isso nao e erro,
  # e so o sinal de cair no caminho de replotagem a partir do CSV.
  runs <- tryCatch(list_runs(base_dir), error = function(e) NULL)
  if (is.null(runs)) return(NULL)
  runs <- runs %>% filter(.data$application | .data$tasks)
  if (nrow(runs) == 0) return(NULL)
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
             cfg = cfg_label(r$runtime, r$scheduler),
             makespan_s = makespan_s)
  })
  dat <- bind_rows(per_run)
  if (nrow(dat) == 0) return(NULL)

  dat <- dat %>%
    mutate(
      # GFLOPS entregues por classe ao longo do run inteiro
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
  write.csv(dat, file.path(plots_dir(), "gpu_split_summary.csv"),
            row.names = FALSE)
  dat
}

#' Le uma decomposicao ja calculada. As colunas do CSV sao exatamente o data
#' frame que as tres figuras consomem, entao replotar dele e equivalente a
#' recalcular -- so nao revalida os rastros.
from_summary <- function(base_dir) {
  f <- file.path(base_dir, "gpu_split_summary.csv")
  if (!file.exists(f)) return(NULL)
  message("sem run dirs; replotando de ", f)
  read.csv(f, stringsAsFactors = FALSE)
}

dat <- from_traces(base_dir)
if (is.null(dat)) dat <- from_summary(base_dir)
if (is.null(dat))
  stop("nem run dirs com parquet nem gpu_split_summary.csv em ", base_dir)

# Notacao e ordem canonicas (plot_style.r): eixo X em "runtime:escalonador" e
# facetas em potrf/getrf_nopiv, casando com o texto do artigo.
dat <- dat %>% mutate(cfg = cfg_factor(.data$cfg), algo = algo_label(.data$algo))
print(as.data.frame(dat), digits = 4)

# CPU/GPU nao sao configuracoes: cor por classe de recurso, com o Set1 pedido
# pelo revisor. O texto do artigo avisa que aqui a cor muda de significado.
cls_fill <- scale_fill_brewer(palette = "Set1")

x_cfg <- theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ---- (1) GFLOPS entregues, empilhado por classe (Figura 3 do artigo) ----
tot <- dat %>% group_by(.data$algo, .data$cfg) %>%
  summarise(total = sum(.data$contrib_gfps), .groups = "drop")
message("para o caption (Fig. contrib): ", base_dir,
        " | total da barra = desempenho do run; ",
        "o segmento de GPU e praticamente igual nas quatro configs")
p1 <- ggplot(dat, aes(.data$cfg, .data$contrib_gfps, fill = .data$class)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f", .data$contrib_gfps)),
            position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 4, colour = "white") +
  geom_text(data = tot, aes(.data$cfg, .data$total,
                            label = sprintf("%.0f", .data$total)),
            inherit.aes = FALSE, vjust = -0.4, size = BASE_SIZE / 3.8,
            fontface = "bold") +
  facet_wrap(~algo) +
  cls_fill +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "GFLOPS entregues") +
  theme_sscad(legend = "bottom") +
  x_cfg
save_plot(p1, "gpu_split_contrib", width = 11, height = 6)

# ---- (2) % dos flops na GPU vs split otimo r* ----
shr <- dat %>% filter(.data$class == "GPU")
message("para o caption (Fig. share): ", base_dir,
        " | losango = r* = cap_GPU/(cap_GPU+cap_CPU) medido no proprio run; ",
        "acima do losango = GPU sobrecarregada")
p2 <- ggplot(shr, aes(.data$cfg, 100 * .data$share_flops)) +
  geom_col(width = 0.6, alpha = 0.9, fill = "#e41a1c") +
  geom_point(aes(y = 100 * .data$rstar), shape = 18, size = 4, colour = "black") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * .data$share_flops)),
            vjust = -0.5, size = BASE_SIZE / 4) +
  facet_wrap(~algo) +
  expand_y_zero() +
  labs(x = NULL, y = "% dos GFLOPS na GPU") +
  theme_sscad(legend = "none") +
  x_cfg
save_plot(p2, "gpu_split_share", width = 11, height = 6)

# ---- (3) capacidade efetiva por classe ----
message("para o caption (Fig. rates): ", base_dir,
        " | GPU = GFLOPS com o device ativo; ",
        "CPU = workers x taxa por worker ocupado")
p3 <- ggplot(dat, aes(.data$cfg, .data$capacity_gfps, fill = .data$class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f", .data$capacity_gfps)),
            position = position_dodge(width = 0.75), vjust = -0.4,
            size = BASE_SIZE / 4) +
  facet_wrap(~algo) +
  cls_fill +
  expand_y_zero() +
  labs(x = NULL, y = "GFLOPS") +
  theme_sscad(legend = "bottom") +
  x_cfg
save_plot(p3, "gpu_split_rates", width = 11, height = 6)
