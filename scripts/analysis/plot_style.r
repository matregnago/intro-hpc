#!/usr/bin/env Rscript
#
# Estilo visual compartilhado por todas as figuras do artigo. Sourced, never run.
#
# Regra de cor (pedido do revisor): o RUNTIME define a matiz -- StarPU sempre
# vermelho, PaRSEC sempre azul -- e o ESCALONADOR varia a tonalidade e o shape.
# Assim o leitor identifica o runtime pela cor mesmo em preto e branco parcial,
# e o escalonador pelo tom/marcador.
#
# Também centraliza: ordem canônica das configs, nomes dos algoritmos, tamanho
# de fonte, largura de canvas e tema. Qualquer figura nova deve usar daqui, para
# que o tamanho da fonte dos elementos gráficos seja o mesmo em todas.

suppressMessages({
  library(ggplot2)
})

# Tamanho de fonte único para todas as figuras do artigo -- em pt SOBRE O CANVAS
# do ggsave, não sobre a página do artigo.
BASE_SIZE <- 14

# Largura do canvas, e por que ela importa tanto quanto o BASE_SIZE: as figuras
# entram no LaTeX escaladas para \linewidth, então a fonte que o leitor enxerga é
#
#   fonte_renderizada = BASE_SIZE * (largura no LaTeX / largura do canvas)
#
# Fixar só o BASE_SIZE não padroniza nada: um canvas mais largo encolhe a fonte
# na mesma proporção. Toda figura nova deve usar uma destas duas larguras.
FIG_WIDTH_IN      <- 11   # padrão (figuras de linha/barra)
FIG_WIDTH_WIDE_IN <- 12   # painéis de rastro, que precisam de mais eixo X

# Ordem canônica: StarPU primeiro, e dentro de cada runtime o escalonador mais
# simples antes. Mesma ordem usada nas linhas dos painéis espaço-tempo.
CONFIG_ORDER <- c("starpu:dmda", "starpu:dmdas", "parsec:lfq", "parsec:gd")

# StarPU = vermelhos, PaRSEC = azuis. Tons escolhidos das escalas sequenciais
# Reds e Blues do ColorBrewer (e não do divergente RdBu, cujo tom claro do
# vermelho puxa para o laranja): dentro de cada runtime só muda a luminosidade,
# como o texto do artigo descreve -- a matiz é o runtime, o tom é o escalonador.
CONFIG_COLORS <- c(
  "starpu:dmda"  = "#a50f15",  # vermelho escuro
  "starpu:dmdas" = "#fb6a4a",  # vermelho claro
  "parsec:lfq"   = "#6baed6",  # azul claro
  "parsec:gd"    = "#08519c"   # azul escuro
)

# Shape redundante com a cor, para leitura em impressão P&B.
CONFIG_SHAPES <- c(
  "starpu:dmda"  = 16,  # círculo
  "starpu:dmdas" = 17,  # triângulo
  "parsec:lfq"   = 18,  # losango
  "parsec:gd"    = 15   # quadrado
)

#' Rótulo canônico de uma configuração: sempre "runtime:escalonador".
#' É a mesma notação usada no texto do artigo -- não divergir.
cfg_label <- function(runtime, scheduler) paste(runtime, scheduler, sep = ":")

#' Converte um vetor de rótulos de config em factor na ordem canônica, mantendo
#' no fim (na ordem em que aparecem) qualquer config fora da lista.
cfg_factor <- function(x) {
  x <- as.character(x)
  extra <- setdiff(unique(x), CONFIG_ORDER)
  factor(x, levels = c(CONFIG_ORDER, sort(extra)))
}

#' Nome canônico do algoritmo para faceta/legenda: sem prefixo de precisão e em
#' minúsculas, casando com o texto do artigo (=potrf=, =getrf_nopiv=).
algo_label <- function(x) {
  v <- tolower(as.character(x))
  sub("^[sdcz](?=(potrf|getrf|geqrf|gemm|syrk|trsm))", "", v, perl = TRUE)
}

#' Escalas de cor/preenchimento/shape por configuração.
scale_config_colour <- function(name = NULL, ...) {
  scale_colour_manual(name = name, values = CONFIG_COLORS, drop = FALSE, ...)
}
scale_config_fill <- function(name = NULL, ...) {
  scale_fill_manual(name = name, values = CONFIG_COLORS, drop = FALSE, ...)
}
scale_config_shape <- function(name = NULL, ...) {
  scale_shape_manual(name = name, values = CONFIG_SHAPES, drop = FALSE, ...)
}

#' Tema único das figuras. `legend` controla a posição -- o padrão é "top", que é
#' onde a legenda de cores fica em todas as figuras do artigo (pedido do
#' revisor); só use outro valor para "none". O resto é fixo, e junto com
#' FIG_WIDTH_IN é o que garante fontes idênticas entre figuras.
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

#' Ancora o eixo Y no zero, sem cortar o topo dos dados.
#' O revisor pediu Y=0 presente para que as alturas relativas sejam legíveis.
#' `limits = c(0, NA)` mantém o topo livre (compatível com scales="free_y"), e a
#' expansão assimétrica evita a faixa negativa que `expand_limits()` deixaria.
expand_y_zero <- function() {
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05)))
}
