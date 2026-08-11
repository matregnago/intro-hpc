#!/usr/bin/env Rscript
#
# Figura 4 do artigo: espaco-tempo StarVZ de StarPU vs PaRSEC no mesmo problema,
# num unico facet_wrap(~cfg) -- StarPU na coluna da esquerda, PaRSEC na direita,
# um escalonador por linha. norm_st_kernels() alinha nomes e cores dos kernels
# entre os dois runtimes.
#
# Uma linha por worker (a agregacao "nodes" do StarVZ colapsa todos os workers de
# uma classe numa lane so e esconde justamente o que a figura quer mostrar). As 4
# configs de uma maquina compartilham o mesmo layout de lanes, o que e o que
# torna correto usar uma unica escala x/y e uma unica legenda.
#
# Uso:  plot_compare_st.r [base_dir]
#   precisa do application.parquet de cada run (rodar parsec_phase1.sh antes,
#   para os runs de GPU do PaRSEC).
# Saida: <PLOTS_DIR>/st_compare_<algo>.{png,pdf}

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

ST_Y_TEXT    <- 9
ST_WIDTH_IN  <- FIG_WIDTH_WIDE_IN
# altura por lane: 0.09in nao da espaco para o rotulo de cada worker; 0.18-0.25
# infla a figura o bastante para o LaTeX flutua-la para uma pagina propria.
ST_UNIT_IN   <- 0.13
# strip, eixo X e margens de UMA linha da grade de facetas.
ST_CHROME_IN <- 1.10

# aggregation "static" = um retangulo por worker (nao por classe); ABE ligado
# (o limite por area sobre o conjunto heterogeneo CPU+GPU e o que a figura
# discute); idleness desligado por causa dos spans assincronos de GPU.
configure <- function(svz) {
  svz$config$st$outliers           <- FALSE
  svz$config$st$aggregation$active <- TRUE
  svz$config$st$aggregation$method <- "static"
  svz$config$st$aggregation$step   <- 100
  svz$config$st$idleness           <- FALSE
  svz$config$st$makespan           <- TRUE
  svz$config$st$labels             <- "ALL"
  svz$config$st$base_size          <- BASE_SIZE
  svz$config$st$abe$active         <- TRUE
  svz$config$st$abe$label          <- TRUE
  svz$config$st$tasks$active       <- FALSE
  svz
}

#' Fim do makespan de um run, para alinhar o eixo X entre os paineis.
run_end <- function(run_dir) {
  ex <- read_exec(run_dir)
  if (is.null(ex) || !nrow(ex)) return(NA_real_)
  max(ex$End, na.rm = TRUE)
}

#' Le e prepara um run: retangulos agregados, lanes, rotulos, makespan e ABE.
prepare_run <- function(run_dir, cfg) {
  message("lendo: ", run_dir)
  svz <- configure(starvz_read(run_dir, selective = FALSE))
  # lanes ANTES de dropar init: o PaRSEC sempre tem um worker de CPU que so
  # executa geracao de matriz (CPU16 no gd, CPU6 no lfq) e sumiria do eixo.
  lanes  <- svz$Application %>%
    distinct(.data$ResourceId, .data$Position, .data$Height) %>%
    arrange(desc(.data$Position))
  labels <- starvz:::yconf(svz$Application, option = svz$config$st$labels)
  svz <- drop_init_kernels(svz)
  svz <- norm_st_kernels(svz)
  # st_time_aggregation()/hl_global_abe() nao sao exportadas, mas sao o que
  # panel_st() chamaria por baixo; acessa-las direto e o que permite juntar os 4
  # runs num so facet_wrap em vez de colar 4 ggplots com patchwork.
  agg <- starvz:::st_time_aggregation(svz$Application, colors = svz$Colors,
                                      step = svz$config$st$aggregation$step) %>%
    mutate(cfg = cfg)
  tend <- max(svz$Application$End)
  # deoverlap_gpu_durations() e obrigatorio antes do LP: com as duracoes de GPU
  # do PaRSEC como vem do rastro, o ABE sai MAIOR que o makespan e deixa de ser
  # um limite inferior. No StarPU o fator e 1 e nada muda.
  abe <- tryCatch(
    starvz:::hl_global_abe(deoverlap_gpu_durations(svz$Application))$Result[[1]],
    error = function(e) {
      message("  ABE falhou: ", conditionMessage(e)); NA_real_
    })
  list(agg = agg, colors = svz$Colors, lanes = lanes, labels = labels, tend = tend, abe = abe)
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

  # Ordem row-major do facet_wrap(ncol=2): linha i pareia starpu[i] com
  # parsec[i], preenchendo com o que sobrar do lado maior.
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

  # Posicoes/alturas sao identicas nas 4 configs, entao a primeira serve para
  # todas. Rotulos de CPU: mantidos. Rotulos de GPU: renomeados GPU0..GPU3 por
  # posicao, porque o texto literal difere (StarPU usa CUDA0_0..3, PaRSEC
  # CUDA0_2..5 -- reserva as streams 0-1 para transferencias) e o facet_wrap usa
  # uma unica escala Y para todos os paineis.
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

  # Equivalentes ao geom_makespan()/geom_abe() do starvz, roteados por facet via
  # a coluna cfg.
  makespan_df <- tibble(
    cfg = factor(cfg_order, levels = cfg_order),
    x   = vapply(prepped, `[[`, numeric(1), "tend"),
    y   = max(lanes$Position) * 0.5
  ) %>% mutate(label = sprintf("%.0f", .data$x))

  abe_df <- tibble(
    cfg  = factor(cfg_order, levels = cfg_order),
    x    = vapply(prepped, `[[`, numeric(1), "abe"),
    ymin = min(lanes$Position),
    ymax = max(lanes$Position) + min(lanes$Height) / 1.25,
    y    = min(lanes$Position) + (max(lanes$Position) - min(lanes$Position)) / 2
  ) %>% mutate(label = sprintf("ABE: %.0f", .data$x)) %>% filter(!is.na(.data$x))

  n  <- unique(na.omit(sub$n)); b <- unique(na.omit(sub$b))
  message(sprintf("para o caption: %s, n=%s, b=%s",
                  alg, paste(n, collapse = "/"), paste(b, collapse = "/")))

  p <- ggplot(agg_all) +
    geom_rect(aes(fill = .data$Task, xmin = .data$Start, xmax = .data$End,
                  ymin = .data$Position + .data$TaskPosition,
                  ymax = .data$Position + .data$TaskPosition + .data$TaskHeight),
              alpha = 0.5) +
    scale_fill_manual(values = fill_colors) +
    geom_segment(data = abe_df,
                aes(x = .data$x, xend = .data$x, y = .data$ymin, yend = .data$ymax),
                inherit.aes = FALSE, linewidth = 5, alpha = 0.7, colour = "gray") +
    geom_text(data = abe_df, aes(x = .data$x, y = .data$y, label = .data$label),
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
