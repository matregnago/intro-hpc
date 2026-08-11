#!/usr/bin/env Rscript
#
# Estilo visual compartilhado por todas as figuras do artigo. Sourced, never run.
#
# Regra de cor: o RUNTIME define a matiz (StarPU vermelho, PaRSEC azul) e o
# ESCALONADOR varia o tom e o shape -- assim o runtime e legivel mesmo em P&B.

suppressMessages({
  library(ggplot2)
})

# Em pt SOBRE O CANVAS do ggsave, nao sobre a pagina do artigo.
BASE_SIZE <- 14

# A largura importa tanto quanto o BASE_SIZE: as figuras entram no LaTeX
# escaladas para \linewidth, entao a fonte que o leitor enxerga e
#
#   fonte_renderizada = BASE_SIZE * (largura no LaTeX / largura do canvas)
#
# Fixar so o BASE_SIZE nao padroniza nada. Toda figura nova usa uma destas duas.
FIG_WIDTH_IN      <- 11   # padrao (figuras de linha/barra)
FIG_WIDTH_WIDE_IN <- 12   # paineis de rastro, que precisam de mais eixo X

# StarPU primeiro, e dentro de cada runtime o escalonador mais simples antes.
CONFIG_ORDER <- c("starpu:dmda", "starpu:dmdas", "parsec:lfq", "parsec:gd")

# Tons das escalas sequenciais Reds/Blues do ColorBrewer (e nao do divergente
# RdBu, cujo vermelho claro puxa para o laranja).
CONFIG_COLORS <- c(
  "starpu:dmda"  = "#a50f15",
  "starpu:dmdas" = "#fb6a4a",
  "parsec:lfq"   = "#6baed6",
  "parsec:gd"    = "#08519c"
)

# Shape redundante com a cor, para leitura em impressao P&B.
CONFIG_SHAPES <- c(
  "starpu:dmda"  = 16,  # circulo
  "starpu:dmdas" = 17,  # triangulo
  "parsec:lfq"   = 18,  # losango
  "parsec:gd"    = 15   # quadrado
)

#' Rotulo canonico de uma configuracao, na mesma notacao usada no texto.
cfg_label <- function(runtime, scheduler) paste(runtime, scheduler, sep = ":")

#' Factor na ordem canonica, mantendo no fim qualquer config fora da lista.
cfg_factor <- function(x) {
  x <- as.character(x)
  extra <- setdiff(unique(x), CONFIG_ORDER)
  factor(x, levels = c(CONFIG_ORDER, sort(extra)))
}

#' Nome canonico do algoritmo: sem prefixo de precisao e em minusculas, casando
#' com o texto do artigo (=potrf=, =getrf_nopiv=).
algo_label <- function(x) {
  v <- tolower(as.character(x))
  sub("^[sdcz](?=(potrf|getrf|geqrf|gemm|syrk|trsm))", "", v, perl = TRUE)
}

scale_config_colour <- function(name = NULL, ...) {
  scale_colour_manual(name = name, values = CONFIG_COLORS, drop = FALSE, ...)
}
scale_config_shape <- function(name = NULL, ...) {
  scale_shape_manual(name = name, values = CONFIG_SHAPES, drop = FALSE, ...)
}

#' Tema unico das figuras; `legend` controla a posicao, o resto e fixo para
#' garantir fontes identicas entre figuras.
theme_sscad <- function(legend = "top") {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      legend.position  = legend,
      legend.title     = element_blank(),
      strip.text       = element_text(face = "bold", size = BASE_SIZE - 1),
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      axis.title       = element_text(size = BASE_SIZE),
      axis.text        = element_text(size = BASE_SIZE - 2),
      panel.grid.minor = element_blank()
    )
}

#' Ancora o eixo Y no zero sem cortar o topo dos dados: limits = c(0, NA) mantem
#' o topo livre (compativel com scales="free_y") e a expansao assimetrica evita
#' a faixa negativa que expand_limits() deixaria.
expand_y_zero <- function() {
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05)))
}
