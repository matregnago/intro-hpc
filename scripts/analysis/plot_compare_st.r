#!/usr/bin/env Rscript
#
# Side-by-side StarVZ space-time for StarPU vs PaRSEC on the same problem, as a
# single facet_wrap(~cfg) plot: StarPU on the left column, PaRSEC on the
# right, one scheduler per row. norm_st_kernels() aligns the kernel
# names/colours across the two runtimes.
#
# One row per worker (not the "nodes" aggregation, which collapses every
# worker of a class into a single occupancy lane -- tried and reverted, it
# hid too much of the per-worker story this figure exists to tell). The 4
# configs of a machine share the exact same lane layout for CPU (same
# ResourceId/Position/Height), which is what makes one shared x/y scale + one
# legend correct.
#
# Y-axis labels: every worker is labelled (StarVZ's labels="ALL"). The one
# thing that does NOT match verbatim across the 4 configs is the GPU
# ResourceId text: StarPU names its 4 GPU workers CUDA0_0..CUDA0_3, PaRSEC
# CUDA0_2..CUDA0_5 (it reserves streams 0-1 for H2D/D2H transfers, which carry
# no compute task and so never appear in Application). Since facet_wrap uses
# one Y scale for every panel, GPU lanes are relabelled generically as
# GPU0..GPU3 (by position, same in all 4 configs); CPU labels are kept as-is
# since those genuinely match.
#
# Usage:  plot_compare_st.r [base_dir]
#   needs each run's application.parquet (run parsec_phase1.sh for the PaRSEC
#   GPU runs first). Outputs <PLOTS_DIR>/st_compare_<algo>.{png,pdf}.

suppressMessages({
  library(tidyverse); library(starvz)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
source(file.path(script_dir, "plot_style.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"

# Meio-termo: 0.09in/unidade (usado quando so um rotulo por painel aparecia)
# nao da espaco para o rotulo de cada worker; 0.1757-0.25 (11-12pt, testados
# antes) da espaco de sobra mas infla a figura o bastante para o LaTeX
# flutua-la para uma pagina propria. 0.13 fica entre os dois.
ST_Y_TEXT    <- 9
ST_WIDTH_IN  <- FIG_WIDTH_WIDE_IN
ST_UNIT_IN   <- 0.13
# Titulo do painel (strip), eixo X e margens de UMA linha da grade de facetas:
# o que cada linha gasta alem das lanes.
ST_CHROME_IN <- 1.10

# panel_st knobs: aggregation on (agrega no tempo para reduzir o numero de
# retangulos, no metodo "static" = um retangulo por worker, nao por classe),
# ABE off (conjunto heterogeneo CPU+GPU), idleness off (spans assincronos de
# GPU), CPB so quando o run tem Dag.
configure <- function(svz) {
  svz$config$st$outliers           <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$method <- "static"
  svz$config$st$aggregation$step   <- 100
  svz$config$st$idleness           <- FALSE
  svz$config$st$makespan           <- TRUE
  svz$config$st$labels             <- "ALL"
  svz$config$st$base_size          <- BASE_SIZE
  svz$config$st$abe$active         <- FALSE
  svz$config$st$abe$label          <- FALSE
  svz$config$st$tasks$active       <- FALSE
  svz$config$st$cpb <- !is.null(svz$Dag) && nrow(svz$Dag) > 0
  svz
}

#' Fim do makespan de um run, para alinhar o eixo X entre os paineis.
run_end <- function(run_dir) {
  ex <- read_exec(run_dir)
  if (is.null(ex) || !nrow(ex)) return(NA_real_)
  max(ex$End, na.rm = TRUE)
}

#' Le e prepara um run: retangulos agregados (com a coluna cfg ja presa), o
#' conjunto COMPLETO de lanes ANTES de dropar kernels de init (o PaRSEC sempre
#' tem um worker de CPU que so executa geracao de matriz -- CPU16 no gd, CPU6
#' no lfq -- que sumiria do eixo sem essa captura), os rotulos (starvz:::yconf,
#' mesma logica que panel_st() usaria), o makespan e o CPB (ou NA se o run nao
#' tem Dag).
prepare_run <- function(run_dir, cfg) {
  message("lendo: ", run_dir)
  svz <- configure(starvz_read(run_dir, selective = FALSE))
  lanes  <- svz$Application %>%
    distinct(.data$ResourceId, .data$Position, .data$Height) %>%
    arrange(desc(.data$Position))
  labels <- starvz:::yconf(svz$Application, option = svz$config$st$labels)
  svz <- drop_init_kernels(svz)
  svz <- norm_st_kernels(svz)
  # st_time_aggregation()/hl_global_cpb() nao sao exportadas pelo starvz, mas
  # sao exatamente o que panel_st() chamaria por baixo dos panos para este run;
  # acessa-las direto e o que permite combinar os 4 runs num so facet_wrap em
  # vez de 4 ggplots colados com patchwork.
  agg <- starvz:::st_time_aggregation(svz$Application, colors = svz$Colors,
                                      step = svz$config$st$aggregation$step) %>%
    mutate(cfg = cfg)
  tend <- max(svz$Application$End)
  cpb  <- if (isTRUE(svz$config$st$cpb)) starvz:::hl_global_cpb(svz)$CPB else NA_real_
  list(agg = agg, colors = svz$Colors, lanes = lanes, labels = labels, tend = tend, cpb = cpb)
}

runs <- list_runs(base_dir) %>% filter(.data$application)
if (nrow(runs) == 0) stop("no runs with application.parquet under ", base_dir,
                          " (run parsec_phase1.sh for the PaRSEC runs)")

for (alg in sort(unique(runs$algo))) {
  sub <- runs %>% filter(.data$algo == alg) %>% arrange(.data$runtime, .data$scheduler)
  starpu <- sub %>% filter(.data$runtime == "starpu")
  parsec <- sub %>% filter(.data$runtime == "parsec")
  if (nrow(starpu) == 0 || nrow(parsec) == 0) {
    message("skip ", alg, ": need both runtimes (have starpu=", nrow(starpu),
            ", parsec=", nrow(parsec), ")")
    next
  }

  # Eixo X comum: o maior makespan do grupo manda em todos os paineis.
  x_end <- suppressWarnings(max(vapply(sub$path, run_end, numeric(1)),
                                na.rm = TRUE))
  if (!is.finite(x_end)) x_end <- NA_real_

  # Ordem row-major (StarPU na coluna 1, PaRSEC na coluna 2 do facet_wrap
  # ncol=2): linha i pareia starpu[i] com parsec[i], preenchendo com o que
  # sobrar do lado maior.
  nrows <- max(nrow(starpu), nrow(parsec))
  reps <- list()
  for (i in seq_len(nrows)) {
    if (i <= nrow(starpu)) reps[[length(reps) + 1]] <-
      list(path = starpu$path[i], cfg = cfg_label(starpu$runtime[i], starpu$scheduler[i]))
    if (i <= nrow(parsec)) reps[[length(reps) + 1]] <-
      list(path = parsec$path[i], cfg = cfg_label(parsec$runtime[i], parsec$scheduler[i]))
  }
  cfg_order <- vapply(reps, function(r) r$cfg, character(1))
  prepped   <- lapply(reps, function(r) prepare_run(r$path, r$cfg))

  agg_all <- bind_rows(lapply(prepped, `[[`, "agg")) %>%
    mutate(cfg = factor(.data$cfg, levels = cfg_order))

  # Posicoes/alturas completas (todos os workers) sao identicas nas 4 configs
  # de uma maquina, entao a primeira ja serve para as 4. Labels de CPU: mantidos
  # como estao (identicos nas 4 configs). Labels de GPU: renomeados GPU0..GPU3
  # por posicao -- ver nota no topo do arquivo sobre por que o rotulo literal
  # (CUDA0_0..3 vs CUDA0_2..5) nao serve para uma unica escala Y compartilhada.
  lanes  <- prepped[[1]]$lanes
  labels <- prepped[[1]]$labels
  labels <- bind_rows(
    labels %>% filter(.data$Height == 1),
    labels %>% filter(.data$Height == 2) %>% arrange(desc(.data$Position)) %>%
      mutate(ResourceId = sprintf("GPU%d", row_number() - 1))
  ) %>% arrange(desc(.data$Position))
  units  <- sum(lanes$Height)

  fill_colors <- setNames(as.character(prepped[[1]]$colors$Color),
                          as.character(prepped[[1]]$colors$Value))

  # Rotulo de makespan por painel (angulo 90, na borda direita de cada run) --
  # equivalente ao geom_makespan() do starvz, mas roteado por facet via a
  # coluna cfg.
  makespan_df <- tibble(
    cfg = factor(cfg_order, levels = cfg_order),
    x   = vapply(prepped, `[[`, numeric(1), "tend"),
    y   = max(lanes$Position) * 0.5
  ) %>% mutate(label = sprintf("%.0f", .data$x))

  # Banda + rotulo de CPB por painel -- equivalente ao geom_cpb() do starvz,
  # mesma logica.
  cpb_df <- tibble(
    cfg  = factor(cfg_order, levels = cfg_order),
    x    = vapply(prepped, `[[`, numeric(1), "cpb"),
    ymin = min(lanes$Position),
    ymax = max(lanes$Position) + min(lanes$Height) / 1.25,
    y    = min(lanes$Position) + (max(lanes$Position) - min(lanes$Position)) / 2
  ) %>% mutate(label = sprintf("CPB: %.0f", .data$x)) %>% filter(!is.na(.data$x))

  n  <- unique(na.omit(sub$n)); b <- unique(na.omit(sub$b))
  message(sprintf("para o caption: %s, n=%s, b=%s",
                  alg, paste(n, collapse = "/"), paste(b, collapse = "/")))

  p <- ggplot(agg_all) +
    geom_rect(aes(fill = .data$Task, xmin = .data$Start, xmax = .data$End,
                  ymin = .data$Position + .data$TaskPosition,
                  ymax = .data$Position + .data$TaskPosition + .data$TaskHeight),
              alpha = 0.5) +
    scale_fill_manual(values = fill_colors) +
    geom_segment(data = cpb_df,
                aes(x = .data$x, xend = .data$x, y = .data$ymin, yend = .data$ymax),
                inherit.aes = FALSE, linewidth = 5, alpha = 0.7, colour = "gray") +
    geom_text(data = cpb_df, aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE, angle = 90, colour = "black", size = BASE_SIZE / 5) +
    geom_text(data = makespan_df, aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE, angle = 90, size = BASE_SIZE / 4) +
    facet_wrap(~cfg, ncol = 2) +
    scale_x_continuous(expand = c(0.05, 0),
                       labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
    scale_y_continuous(breaks = labels$Position, labels = labels$ResourceId,
                       expand = c(0.05, 0)) +
    coord_cartesian(xlim = c(0, x_end), ylim = c(0, NA)) +
    labs(x = "Time [ms]", y = "Application Workers") +
    theme_sscad() +
    theme(axis.text.y = element_text(size = ST_Y_TEXT), panel.grid = element_blank())

  save_plot(p, paste0("st_compare_", alg), width = ST_WIDTH_IN,
            height = (units * ST_UNIT_IN + ST_CHROME_IN) * nrows)
}
