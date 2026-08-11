#!/usr/bin/env Rscript
#
# Figura 3 do artigo: para cada run de GPU, decompoe o trabalho da fatoracao por
# classe de recurso (CPU vs GPU) e mostra quanto cada uma entregou de GFLOPS
# (gflop da classe / makespan). A soma dos segmentos e o desempenho do run.
#
# Flops por tarefa vem da mesma formula nos dois runtimes (kernel x b do nome do
# run), e nao da coluna GFlop do StarVZ, que o PaRSEC nao tem.
#
# Uso:  plot_gpu_split.r [poti_dir tupi_dir]
#   com base_dirs -> recalcula dos rastros e reescreve o summary.
#   sem args e com plots/final/gpu_split_summary.csv -> so replota do CSV.
#   default: data/traces_poti_*/runs e data/traces_tupi_*/runs mais recentes.
#
# Saida: <PLOTS_DIR>/gpu_split_contrib.{png,pdf} + gpu_split_summary.csv

suppressMessages({
  library(arrow); library(dplyr); library(ggplot2)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
source(file.path(script_dir, "plot_style.r"))

args <- commandArgs(trailingOnly = TRUE)
dirs <- if (length(args) >= 2) {
  args[1:2]
} else if (length(args) == 1) {
  stop("passe 0 ou 2 base_dirs (poti e tupi); recebeu 1")
} else {
  poti_cands <- sort(Sys.glob("data/traces_poti_*/runs"), decreasing = TRUE)
  tupi_cands <- sort(Sys.glob("data/traces_tupi_*/runs"), decreasing = TRUE)
  if (length(poti_cands) == 0) stop("nenhum data/traces_poti_*/runs encontrado")
  if (length(tupi_cands) == 0) stop("nenhum data/traces_tupi_*/runs encontrado")
  c(poti_cands[1], tupi_cands[1])
}

# flops em unidades de b^3
FLOP_COEF <- c(gemm = 2, trsm = 1, syrk = 1, herk = 1, trmm = 1, potrf = 1 / 3,
               getrf_nopiv = 2 / 3, getrf = 2 / 3)

node_of <- function(base_dir) {
  m <- regexpr("poti|tupi", base_dir)
  if (m < 0) stop("nao foi possivel identificar poti/tupi em ", base_dir)
  regmatches(base_dir, m)
}
gpu_label <- function(node) switch(node,
  poti = "poti (RTX 4070)", tupi = "tupi (RTX 4090)", node)

#' Recalcula a decomposicao de UM no a partir dos rastros de cada run.
from_traces <- function(base_dir) {
  # list_runs() aborta quando nao ha run dirs; aqui isso nao e erro, e so o
  # sinal de cair no caminho de replotagem a partir do CSV.
  runs <- tryCatch(list_runs(base_dir), error = function(e) NULL)
  if (is.null(runs)) return(NULL)
  runs <- runs %>% filter(.data$application | .data$tasks)
  if (nrow(runs) == 0) return(NULL)
  message("runs: ", paste(runs$dir, collapse = ", "))

  per_run <- lapply(seq_len(nrow(runs)), function(i) {
    r  <- runs[i, ]
    ex <- read_exec(r$path)
    if (is.null(ex)) return(NULL)
    # kernels de init rodam antes da regiao cronometrada
    ex <- ex %>%
      filter(!.data$kernel %in% INIT_KERNELS,
             .data$kernel %in% names(FLOP_COEF)) %>%
      mutate(class = ifelse(grepl("CUDA", .data$worker), "GPU", "CPU"),
             gflop = FLOP_COEF[.data$kernel] * as.numeric(r$b)^3 / 1e9)
    if (nrow(ex) == 0) return(NULL)
    makespan_s <- (max(ex$End) - min(ex$Start)) / 1000
    ex %>%
      group_by(.data$class) %>%
      # active_s pela uniao dos intervalos: as streams de GPU do PaRSEC se
      # sobrepoem, entao somar Duration inflaria o tempo ativo do device.
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

  node    <- node_of(base_dir)
  gpu_lbl <- gpu_label(node)
  dat %>% mutate(node = node, gpu = gpu_lbl,
                 contrib_gfps = .data$gflop / .data$makespan_s)
}

#' Le uma decomposicao ja calculada (as duas maquinas juntas): as colunas do CSV
#' sao o proprio data frame que a figura consome.
from_summary <- function() {
  f <- file.path(plots_dir(), "gpu_split_summary.csv")
  if (!file.exists(f)) return(NULL)
  message("sem run dirs; replotando de ", f)
  read.csv(f, stringsAsFactors = FALSE)
}

dat <- bind_rows(lapply(dirs, from_traces))
if (nrow(dat) == 0) dat <- from_summary()
if (is.null(dat) || nrow(dat) == 0)
  stop("nem run dirs com parquet em ", paste(dirs, collapse = " e "),
      " nem gpu_split_summary.csv em ", plots_dir())
write.csv(dat, file.path(plots_dir(), "gpu_split_summary.csv"), row.names = FALSE)

dat <- dat %>% mutate(cfg  = cfg_factor(.data$cfg),
                      algo = algo_label(.data$algo),
                      node = factor(.data$node, levels = c("poti", "tupi")))
print(as.data.frame(dat), digits = 4)

tot <- dat %>% group_by(.data$algo, .data$gpu, .data$cfg) %>%
  summarise(total = sum(.data$contrib_gfps), .groups = "drop")
message("para o caption: ", paste(dirs, collapse = ", "),
        " | total da barra = desempenho do run")

# CPU/GPU nao sao configuracoes: aqui a cor e a classe de recurso, nao a config.
# scales="free" porque os tetos de FP64 da 4070 e da 4090 diferem demais.
p <- ggplot(dat, aes(.data$cfg, .data$contrib_gfps, fill = .data$class)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f", .data$contrib_gfps)),
            position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 4, colour = "white") +
  geom_text(data = tot, aes(.data$cfg, .data$total,
                            label = sprintf("%.0f", .data$total)),
            inherit.aes = FALSE, vjust = -0.4, size = BASE_SIZE / 3.8,
            fontface = "bold") +
  facet_grid(algo ~ gpu, scales = "free") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "GFLOPS entregues") +
  theme_sscad() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p, "gpu_split_contrib", width = FIG_WIDTH_IN, height = 9)
