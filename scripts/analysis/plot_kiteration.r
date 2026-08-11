#!/usr/bin/env Rscript
#
# Figura 5 do artigo: k-iteration (iteracoes no Y, tempo no X) via
# panel_kiteration() do StarVZ, para comparar a "assinatura" dos escalonadores ao
# longo do laco externo da fatoracao.
#
# O PaRSEC nao tem Application$Iteration no rastro; parsec_tasks_to_parquet.r a
# deriva. Como panel_kiteration so precisa de $Application/$Colors/$config, uma
# lista starvz_data montada na mao basta (sem starvz_read).
#
# Uso:  plot_kiteration.r <run_dir|base_dir> [run_dir ...]
#   run dirs precisam de tasks.parquet com coluna Iteration nao-NA.
# Saida: uma figura por algoritmo, todas as configs no mesmo eixo de tempo:
#   <PLOTS_DIR>/kiteration_compare_<algo>.{png,pdf}

suppressMessages({
  library(ggplot2); library(starvz); library(patchwork)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))
source(file.path(script_dir, "plot_style.r"))

# Tabela Colors no formato que panel_kiteration espera (Value character, sem Use).
build_colors <- function(values) {
  vals <- sort(unique(values))
  cols <- KERNEL_COLORS[vals]
  missing <- is.na(cols)
  if (any(missing)) cols[missing] <- scales::hue_pal()(sum(missing))
  tibble(Value = vals, Color = unname(cols))
}

# panel_kiteration de um run; NULL se o run nao tem dado utilizavel.
build_panel <- function(run_dir) {
  meta <- parse_run_id(run_dir)
  tasks <- read_tasks_norm(run_dir)            # StartTime/EndTime -> ms
  if (is.null(tasks)) { message("no tasks.parquet in ", run_dir, ", skip"); return(NULL) }
  if (!"Iteration" %in% names(tasks) || all(is.na(tasks$Iteration))) {
    message("no usable Iteration in ", run_dir,
            " (re-run parsec_tasks_to_parquet.r), skip"); return(NULL)
  }

  app <- tasks %>%
    filter(!is.na(.data$StartTime), !is.na(.data$EndTime), !is.na(.data$Iteration)) %>%
    transmute(
      Iteration    = as.integer(.data$Iteration),
      Subiteration = 0L,
      Node         = 0L,
      Start        = .data$StartTime,
      End          = .data$EndTime,
      Value        = norm_kernel(.data$Name)
    ) %>%
    # so trabalho de fatoracao: names(KERNEL_COLORS) e exatamente esse conjunto
    # apos o alias (geracao de matriz e sinks do runtime caem fora).
    filter(.data$Value %in% names(KERNEL_COLORS))
  if (nrow(app) == 0) { message("no timed iterated tasks in ", run_dir, ", skip"); return(NULL) }

  svz <- structure(list(
    Application = app,
    Colors      = build_colors(app$Value),
    config = list(
      base_size = BASE_SIZE, expand = 0.05,
      limits = list(start = NA, end = NA),
      # subite/pernode precisam ser FALSE explicito: panel_kiteration erra em NULL.
      kiteration = list(subite = FALSE, pernode = FALSE,
                        middlelines = NULL, legend = TRUE)
    )
  ), class = "starvz_data")

  label <- meta$label
  p <- panel_kiteration(svz) + ggtitle(label) + ylab("Iteration (k)")
  message(sprintf("    %s [%s]: %d tasks over %d iterations",
                  label, meta$algo, nrow(app), max(app$Iteration) + 1L))
  list(panel = p, label = label,
       runtime   = ifelse(is.na(meta$runtime), "run", meta$runtime),
       scheduler = ifelse(is.na(meta$scheduler), basename(run_dir), meta$scheduler),
       algo      = ifelse(is.na(meta$algo), "run", meta$algo),
       x_max = max(app$End))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: plot_kiteration.r <run_dir|base_dir> [run_dir ...]")

# Cada arg e um run dir (tem tasks.parquet) ou um base dir a expandir.
expand_arg <- function(a) {
  if (file.exists(file.path(a, "tasks.parquet"))) return(a)
  subs <- list.dirs(a, recursive = FALSE)
  subs[file.exists(file.path(subs, "tasks.parquet"))]
}
run_dirs <- unique(unlist(lapply(args, expand_arg)))
if (length(run_dirs) == 0)
  stop("no run dirs with tasks.parquet under: ", paste(args, collapse = ", "))

panels <- Filter(Negate(is.null), lapply(run_dirs, build_panel))
if (length(panels) == 0) stop("no runs with usable k-iteration data")

config_rank <- function(p) {
  m <- match(cfg_label(p$runtime, p$scheduler), CONFIG_ORDER)
  if (is.na(m)) length(CONFIG_ORDER) + 1L else m
}

# Uma figura por algoritmo: dentro de um algo todos os runs compartilham N/b,
# entao um eixo de tempo comum mantem os makespans comparaveis.
by_algo <- split(panels, vapply(panels, function(p) p$algo, character(1)))
for (algo in names(by_algo)) {
  grp <- by_algo[[algo]]
  grp <- grp[order(vapply(grp, config_rank, integer(1)))]
  x_end <- max(vapply(grp, function(p) p$x_max, numeric(1)))
  # Legenda no topo, como nas demais figuras. Nao basta legend.position="top": o
  # patchwork encaixa o guia coletado ABAIXO da linha dos ggtitle(), em cima do
  # titulo do painel da direita. guide_area() reserva uma faixa propria acima de
  # tudo, e heights da a ela so o que precisa.
  combined <- (guide_area() / wrap_plots(
    lapply(grp, function(p) p$panel + coord_cartesian(xlim = c(0, x_end))),
    ncol = 2
  )) +
    plot_layout(guides = "collect", heights = c(1, 22)) &
    theme(legend.position = "top")
  # Canvas em FIG_WIDTH_WIDE_IN (era 16x4.5 por linha): mesma razao de aspecto,
  # mesma area na pagina, mas a fonte renderizada sobe ~33% -- este era o painel
  # com a menor fonte do artigo.
  save_plot(combined, paste0("kiteration_compare_", algo),
            width = FIG_WIDTH_WIDE_IN, height = 3.375 * ceiling(length(grp) / 2))
}
