library(ggplot2)
library(dplyr)
library(readr)

results <- read_csv("data/chameleon-runtimes_782024/results.csv",
                    show_col_types = FALSE)

agg <- results |>
  group_by(runtime, scheduler, algorithm, n) |>
  summarise(gflops = mean(gflops), .groups = "drop") |>
  mutate(config = paste(runtime, scheduler, sep = "/"))

p <- ggplot(agg, aes(x = n, y = gflops, color = config)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ algorithm) +
  labs(x = "N (matriz N x N)", y = "GFLOPS", color = NULL,
       title = "Chameleon CPU-only (cei) - b=480, t=24") +
  theme_bw()

dir.create("plots", showWarnings = FALSE)
ggsave("plots/gflops_vs_n_r.png", p, width = 10, height = 4, dpi = 140)
ggsave("plots/gflops_vs_n_r.pdf", p, width = 10, height = 4)
