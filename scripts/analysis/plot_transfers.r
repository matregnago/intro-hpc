#!/usr/bin/env Rscript
#
# Painel de movimentacao de dados H2D/D2H a partir dos contadores nativos que
# scripts/run.sh (TRACE_STATS=1) grava no .log de cada run -- funciona nos dois
# runtimes, sem depender do trace:
#   StarPU: blocos "Data transfer stats" (STARPU_BUS_STATS) + "Worker stats".
#   PaRSEC: linha "|All Devs |" (device_show_statistics).
#
# Uso:  plot_transfers.r [base_dir] [algo] [maquina]
#   algo (opcional) filtra os runs; maquina rotula o painel (default extraido
#   do base_dir). Saida: plots/transfers_d2h.{png,pdf} + transfers_summary.csv.

suppressMessages({
  library(dplyr); library(ggplot2); library(stringr); library(tidyr)
})

this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "."
source(file.path(script_dir, "trace_common.r"))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "data/manual_traces_20260621_172820/runs"
algo_sel <- if (length(args) >= 2) args[[2]] else NA
machine  <- if (length(args) >= 3) args[[3]] else
  str_to_title(str_extract(base_dir, "poti|tupi"))

num <- function(x) as.numeric(str_replace_all(x, "[^0-9.eE+-]", ""))

# StarPU: "Data transfer stats" + "Worker stats" (GB -> MB).
parse_starpu <- function(lines) {
  h2d <- str_match(lines, "NUMA \\d+ -> CUDA \\d+\\s+([0-9.eE+-]+) GB.*transfers : (\\d+)")
  d2h <- str_match(lines, "CUDA \\d+ -> NUMA \\d+\\s+([0-9.eE+-]+) GB.*transfers : (\\d+)")
  h2d <- h2d[!is.na(h2d[, 1]), , drop = FALSE]
  d2h <- d2h[!is.na(d2h[, 1]), , drop = FALSE]
  if (nrow(h2d) == 0 && nrow(d2h) == 0) return(NULL)
  # CUDA tasks: the "N task(s)" line(s) following a "CUDA <dev>" worker line.
  cuda_tasks <- 0L
  w <- grep("Worker stats", lines)
  if (length(w)) {
    blk <- lines[w[1]:length(lines)]
    idx <- grep("^CUDA ", str_trim(blk))
    for (i in idx) {
      m <- str_match(blk[i + 1], "(\\d+) task")
      if (!is.na(m[1])) cuda_tasks <- cuda_tasks + as.integer(m[2])
    }
  }
  tibble(
    h2d_mb      = sum(num(h2d[, 2]) * 1000),
    d2h_mb      = sum(num(d2h[, 2]) * 1000),
    d2h_xfers   = sum(as.integer(d2h[, 3])),
    gpu_kernels = cuda_tasks
  )
}

# PaRSEC: "|All Devs |" summary. PaRSEC scales the unit with the volume
# ("21.09MB(33.33)", "61.40GB(100.00)"), so the parse must be unit-aware.
parse_parsec <- function(lines) {
  ln <- grep("\\|All Devs", lines, value = TRUE)
  if (!length(ln)) return(NULL)
  f <- str_split(ln[1], "\\|")[[1]] |> str_trim()
  # f: "", "All Devs", kernels, %, H2D_req, H2D_xfer(%), D2H_req, D2H_xfer(%), ""
  mb <- function(s) {
    m <- str_match(s, "([0-9.]+)\\s*([KMGT]?)B")
    mult <- c("1e-6", K = "1e-3", M = "1", G = "1e3", T = "1e6")
    key  <- ifelse(m[, 3] == "", 1L, match(m[, 3], names(mult)))
    as.numeric(m[, 2]) * as.numeric(mult[key])
  }
  tibble(
    h2d_mb      = mb(f[6]),
    d2h_mb      = mb(f[8]),
    d2h_xfers   = NA_real_,                       # PaRSEC reports MB, not op count
    gpu_kernels = as.integer(num(f[3]))
  )
}

runs <- list_runs(base_dir)
if (!is.na(algo_sel)) runs <- runs %>% filter(.data$algo == algo_sel)
rows <- lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  log <- list.files(r$path, pattern = "\\.log$", full.names = TRUE)
  if (!length(log)) return(NULL)
  lines <- readLines(log[1], warn = FALSE)
  st <- if (r$runtime == "starpu") parse_starpu(lines) else parse_parsec(lines)
  if (is.null(st)) { message("  sem stats em ", r$dir); return(NULL) }
  bind_cols(tibble(algo = r$algo, runtime = r$runtime,
                   cfg = paste0(r$runtime, ":", r$scheduler)), st)
})
res <- bind_rows(rows)
if (nrow(res) == 0) stop("nenhum log com contadores de transferencia em ", base_dir,
                         " -- rodar a coleta com TRACE_STATS=1")

res <- res %>% mutate(d2h_mb_per_kernel = .data$d2h_mb / .data$gpu_kernels)
print(res, n = Inf)
write.csv(res, file.path(plots_dir(), "transfers_summary.csv"), row.names = FALSE)

long <- res %>%
  pivot_longer(c("h2d_mb", "d2h_mb"), names_to = "dir", values_to = "mb") %>%
  mutate(dir = factor(ifelse(.data$dir == "h2d_mb", "H2D", "D2H"),
                      levels = c("H2D", "D2H")))

# com um unico algoritmo a faixa do facet vira o nome da maquina
one_algo <- n_distinct(long$algo) == 1
long <- long %>%
  mutate(panel = if (one_algo && !is.na(machine)) machine else .data$algo)

p <- ggplot(long, aes(.data$cfg, .data$mb, fill = .data$dir)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f", .data$mb / 1000)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  facet_wrap(~panel) +
  scale_y_continuous(labels = function(x) sprintf("%.0f", x / 1000)) +
  labs(title = "Trafego PCIe H2D/D2H por runtime:scheduler",
       x = NULL, y = "volume transferido (GB)", fill = "direcao") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
save_plot(p, "transfers_d2h", width = 11, height = 6)
