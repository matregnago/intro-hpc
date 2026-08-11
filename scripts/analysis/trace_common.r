#!/usr/bin/env Rscript
#
# Helpers compartilhados pelas figuras de comparacao entre runtimes: parse dos
# nomes de run dir, descoberta de runs, normalizacao dos nomes de kernel entre
# StarPU/PaRSEC e escrita das figuras. Sourced, never run.

suppressMessages({
  library(dplyr); library(stringr); library(tibble); library(arrow)
})

# Hue estavel por kernel, para que a mesma operacao mantenha a cor entre runs,
# runtimes e escalonadores.
KERNEL_COLORS <- c(
  potrf = "#e41a1c", trsm = "#377eb8", syrk = "#4daf4a", gemm = "#984ea3",
  herk = "#4daf4a", trmm = "#377eb8",
  geqrt = "#e41a1c", ormqr = "#984ea3", qr_couple = "#377eb8", qr_apply = "#4daf4a",
  getrf = "#e41a1c", getrf_nopiv = "#e41a1c"
)

#' Tabela Colors no estilo StarVZ: hue de KERNEL_COLORS, fallback deterministico
#' via hue_pal() para kernels desconhecidos.
KERNEL_COLORS_TO_TIBBLE <- function(values) {
  vals  <- sort(unique(as.character(values)))
  cols  <- unname(KERNEL_COLORS[norm_kernel(vals)])
  miss  <- is.na(cols)
  if (any(miss)) cols[miss] <- scales::hue_pal()(sum(miss))
  tibble(Value = factor(vals, levels = vals), Color = cols, Use = TRUE)
}

#' Remove as tarefas de geracao/init da matriz. Os rastros FxT do StarPU nao as
#' contem, entao num painel lado a lado elas so apareceriam na coluna do PaRSEC.
#' Chamar ANTES de norm_st_kernels().
drop_init_kernels <- function(svz) {
  if (is.null(svz$Application) || !"Value" %in% names(svz$Application)) return(svz)
  keep <- !norm_kernel(svz$Application$Value) %in% INIT_KERNELS
  svz$Application <- svz$Application[keep, ]
  svz
}

#' Normaliza Application$Value e reconstroi Colors, para que o painel pinte o
#' mesmo kernel com a mesma cor no StarPU ("dgemm") e no PaRSEC ("gemm").
norm_st_kernels <- function(svz) {
  if (is.null(svz$Application) || !"Value" %in% names(svz$Application)) return(svz)
  app <- svz$Application
  app$Value <- norm_kernel(app$Value)
  if ("lowercase" %in% names(app))
    app$lowercase <- tolower(as.character(app$Value))
  svz$Application <- app
  svz$Colors <- KERNEL_COLORS_TO_TIBBLE(app$Value)
  svz
}

# Kernels de computacao (sem prefixo de precisao). Os codelets de apply do QR do
# Chameleon terminam em 't': o StarPU emite o esquema TS/TT e o PaRSEC o TP.
CHOLESKY_KERNELS <- c("potrf", "trsm", "syrk", "gemm", "herk", "trmm")
QR_KERNELS       <- c("geqrt", "ormqr", "unmqr",
                      "tsqrt", "ttqrt", "tpqrt",        # couple factor
                      "tsmqrt", "ttmqrt", "tpmqrt")     # couple apply
LU_KERNELS       <- c("getrf", "getrf_nopiv")
INIT_KERNELS     <- c("plgsy", "plrnt", "plghe", "lacpy", "laset")
COMPUTE_KERNELS  <- c(CHOLESKY_KERNELS, QR_KERNELS, LU_KERNELS, INIT_KERNELS,
                      "qr_couple", "qr_apply")

# Mesma matematica sob nomes diferentes: ormqr = unmqr, e os esquemas TS/TT do
# StarPU e TP do PaRSEC colapsam num rotulo por estagio do QR.
KERNEL_ALIAS <- c(
  unmqr  = "ormqr",
  tsqrt  = "qr_couple",  ttqrt  = "qr_couple",  tpqrt  = "qr_couple",
  tsmqrt = "qr_apply",   ttmqrt = "qr_apply",   tpmqrt = "qr_apply"
)

#' Colapsa o nome de um kernel num rotulo comum: tira o prefixo de precisao
#' d/s/c/z ("dgemm" -> "gemm") e aplica KERNEL_ALIAS.
norm_kernel <- function(value) {
  v <- tolower(as.character(value))
  stem <- ifelse(
    nchar(v) > 1L &
      substr(v, 1, 1) %in% c("d", "s", "c", "z") &
      substr(v, 2, nchar(v)) %in% COMPUTE_KERNELS,
    substr(v, 2, nchar(v)), v
  )
  ifelse(stem %in% names(KERNEL_ALIAS), KERNEL_ALIAS[stem], stem)
}

#' Parse do basename de um run dir
#' ("<idx>_<runtime>_<sched>_<algo>_n<N>_b<B>_rep<R>", o esquema que
#' scripts/run.sh escreve); campos que nao casam voltam NA.
parse_run_id <- function(dir_name) {
  d <- basename(dir_name)
  m <- str_match(
    d,
    "^(\\d+)_([a-z]+)_([a-z]+)_([a-z][a-z0-9_]*)_n(\\d+)_b(\\d+)_rep(\\d+)"
  )
  tibble(
    dir       = d,
    idx       = m[, 2],
    runtime   = m[, 3],
    scheduler = m[, 4],
    algo      = m[, 5],
    n         = as.integer(m[, 6]),
    b         = as.integer(m[, 7]),
    rep       = as.integer(m[, 8]),
    label     = ifelse(is.na(m[, 3]), d,
                       paste0(m[, 3], ":", m[, 4]))
  )
}

#' Quais dos parquets relevantes um run dir tem.
run_has <- function(run_dir) {
  has <- function(f) file.exists(file.path(run_dir, f))
  list(
    application = has("application.parquet"),
    tasks       = has("tasks.parquet"),
    dag         = has("dag.parquet"),
    variable    = has("variable.parquet"),
    states      = has("states.parquet")
  )
}

#' Descobre os runs sob um base dir (o runs/ de um job): uma linha por run, com
#' as colunas de parse_run_id() + `path` + os flags de run_has().
list_runs <- function(base_dir, pattern = NULL) {
  dirs <- list.dirs(base_dir, recursive = FALSE)
  if (!is.null(pattern)) dirs <- dirs[str_detect(basename(dirs), pattern)]
  if (length(dirs) == 0) {
    stop("no run dirs under ", base_dir,
         if (!is.null(pattern)) paste0(" matching /", pattern, "/") else "")
  }
  meta <- bind_rows(lapply(dirs, parse_run_id))
  meta$path <- dirs
  avail <- bind_rows(lapply(dirs, function(d) as_tibble(run_has(d))))
  bind_cols(meta, avail) %>% arrange(.data$dir)
}

#' Tamanho da uniao dos intervalos [start, end]. Imune a sobreposicao: os spans
#' de GPU do PaRSEC sao observados no host e se sobrepoem entre streams, entao
#' uma soma simples inflaria o tempo ativo do device.
union_length <- function(start, end) {
  stopifnot(length(start) == length(end))
  keep <- !is.na(start) & !is.na(end) & end > start
  s <- start[keep]; e <- end[keep]
  if (length(s) == 0) return(0)
  o   <- order(s)
  s   <- s[o]; e <- e[o]
  cme <- cummax(e)                        # max end ate i
  brk <- c(TRUE, s[-1] > cme[-length(s)]) # TRUE onde abre um novo componente
  idx <- which(brk)
  fin <- c(idx[-1] - 1, length(s))
  sum(cme[fin] - s[idx])
}

#' Desinfla as duracoes de GPU auto-sobrepostas de uma tabela Application.
#'
#' Como os spans do PaRSEC sao cronometrados no HOST (submit -> sync), varias
#' tarefas da MESMA stream parecem rodar ao mesmo tempo: sum(Duration) de uma
#' lane passa de 2.6-3.5x o tempo que ela ficou de fato ocupada. Dividir cada
#' tarefa pelo fator sum/union da lane transforma cada span na sua fatia
#' serial-equivalente, que e o que uma metrica por area (ABE) precisa para
#' continuar sendo um limite inferior. No StarPU o fator e 1 e isso e um no-op.
deoverlap_gpu_durations <- function(app) {
  if (is.null(app) || !all(c("ResourceId", "Start", "End", "Duration") %in% names(app)))
    return(app)
  gpu <- grepl("CUDA", app$ResourceId)
  if (!any(gpu)) return(app)
  fac <- app[gpu, ] %>%
    group_by(.data$ResourceId) %>%
    summarise(f = sum(.data$Duration) / union_length(.data$Start, .data$End),
              .groups = "drop")
  f <- setNames(fac$f, fac$ResourceId)[as.character(app$ResourceId)]
  f[is.na(f) | !is.finite(f) | f < 1] <- 1
  app$Duration <- app$Duration / f
  app
}

#' Escala que leva os tempos de um run para MILISSEGUNDOS. Os conversores
#' offline do PaRSEC -- os unicos que escrevem states.parquet -- emitem
#' microssegundos; os parquets da fase 1 do StarVZ ja estao em ms.
time_scale_to_ms <- function(run_dir) {
  if (file.exists(file.path(run_dir, "states.parquet"))) 1e-3 else 1
}

#' tasks.parquet de um run com StartTime/EndTime normalizados para ms.
read_tasks_norm <- function(run_dir) {
  f <- file.path(run_dir, "tasks.parquet")
  if (!file.exists(f)) return(NULL)
  s <- time_scale_to_ms(run_dir)
  read_parquet(f) %>%
    mutate(StartTime = .data$StartTime * s, EndTime = .data$EndTime * s)
}

#' Intervalos de execucao por kernel de um run, agnostico de runtime: prefere
#' application.parquet e cai para tasks.parquet (a fonte que e completa para a
#' GPU do PaRSEC). Tempos em ms.
#' @return tibble(kernel, Start, End, Duration, worker, src) so de
#'   COMPUTE_KERNELS, ou NULL se nenhuma fonte tem estados utilizaveis.
read_exec <- function(run_dir) {
  pick <- function(tb) tb %>%
    filter(.data$kernel %in% COMPUTE_KERNELS, !is.na(.data$Duration),
           .data$Duration > 0)

  appf <- file.path(run_dir, "application.parquet")
  if (file.exists(appf)) {
    a <- read_parquet(appf)
    if (nrow(a) > 0 && "Value" %in% names(a)) {
      out <- pick(tibble(kernel = norm_kernel(a$Value), Start = a$Start,
                         End = a$End, Duration = a$Duration,
                         worker = as.character(a$ResourceId), src = "application"))
      if (nrow(out) > 0) return(out)
    }
  }
  tkf <- file.path(run_dir, "tasks.parquet")
  if (file.exists(tkf)) {
    s <- time_scale_to_ms(run_dir)
    t <- read_parquet(tkf)
    out <- pick(tibble(kernel = norm_kernel(t$Name), Start = t$StartTime * s,
                       End = t$EndTime * s, Duration = (t$EndTime - t$StartTime) * s,
                       worker = as.character(t$WorkerId), src = "tasks"))
    if (nrow(out) > 0) return(out)
  }
  NULL
}

#' Diretorio de saida das figuras/tabelas; honra $PLOTS_DIR (default "plots").
plots_dir <- function() {
  d <- Sys.getenv("PLOTS_DIR", "plots")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

#' Salva PNG (300 dpi) + PDF (cairo) sob plots_dir().
save_plot <- function(plot, stem, width = 12, height = 6) {
  od  <- plots_dir()
  png <- file.path(od, paste0(stem, ".png"))
  pdf <- file.path(od, paste0(stem, ".pdf"))
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 300,
                  limitsize = FALSE)
  ggplot2::ggsave(pdf, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, limitsize = FALSE)
  message("salvo: ", png, " e ", pdf)
}
