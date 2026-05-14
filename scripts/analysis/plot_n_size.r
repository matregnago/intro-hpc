library(ggplot2)
library(dplyr)
library(readr)

results <- read_csv("data/chameleon-runtimes_782087/results.csv",
  show_col_types = FALSE
)

agg <- results |>
  group_by(runtime, scheduler, algorithm, n) |>
  summarise(
    gflops_mean = mean(gflops),
    gflops_min = min(gflops),
    gflops_max = max(gflops),
    .groups = "drop"
  ) |>
  mutate(config = paste(runtime, scheduler, sep = "/"))

p <- ggplot(agg, aes(x = n, y = gflops_mean, color = config)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(
    ymin = gflops_min,
    ymax = gflops_max
  ), width = 0.02 * max(agg$n)) +
  facet_wrap(~algorithm) +
  labs(
    x = "N (matriz N x N)", y = "GFLOPS", color = NULL,
    title = "Chameleon CPU-only (cei) - b=480, t=24"
  ) +
  theme_bw()

dir.create("plots", showWarnings = FALSE)
ggsave("plots/gflops_vs_n_r.png", p, width = 10, height = 4, dpi = 140)
ggsave("plots/gflops_vs_n_r.pdf", p, width = 10, height = 4)
