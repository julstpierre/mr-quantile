#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(TwoSampleMR)

# Load custom estimator functions
source("utils.R")

# Load data
# D <- sim_indv_dat(list(p = 100, theta = 0, n = 50000, q = 0.6))
chunk <- D_S32[1:400]

# Compare methods
tag_method <- function(res_list, method) {
  bind_rows(res_list) %>% 
    mutate(method = method)
}

res_all <- bind_rows(
  tag_method(lapply(chunk, fit_mr_median), "mr_median"),
  tag_method(lapply(chunk, fit_mr_quantreg), "mr_quantreg"),
  tag_method(lapply(chunk, fit_mr_quantreg_ma, weights = "bic"), "mr_quantreg_ma")
)

# Summarize
res_all %>%
  group_by(method) %>%
  mutate(theta = chunk[[1]]$theta) %>%
  summarise(
    bias = mean(Estimate - theta),
    sd_estimate = sd(Estimate),
    mean_se = mean(StdError),
    rmse = sqrt(mean((Estimate - theta)^2)),
    coverage = mean((LCL <= theta) & (UCL >= theta)),
    power = mean((LCL > 0) | (UCL < 0))
  )

# Plot bias
mutate(res_all, theta = chunk[[1]]$theta) %>%
  mutate(bias = Estimate - theta) %>%
  ggplot(aes(x = method, y = bias, fill = method)) +
  geom_boxplot() +
  # coord_cartesian(ylim = c(-0.05, 0.05)) +
  theme_minimal(base_size = 13) +
  labs(
    title = bquote("Bias distribution, n = " * .(50000) *
                     ", p = " * .(50) * ", Invalid IVs = " * .(0.6 * 100) * "%, " *
                     theta == .(0)
    ),
    y = "Bias",
    x = "Method"
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# -----------------------------
# Save results
# -----------------------------
