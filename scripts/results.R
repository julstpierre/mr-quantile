#!/usr/bin/env Rscript

# -----------------------------
# Load packages
# -----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(ggplot2)

# -----------------------------
# Parameters
# -----------------------------
results_dir <- "../results/rds"
plot_dir <- "../results/plots"
scenario_prefix <- "D_S"

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# Scenario metadata
# -----------------------------
scenario_metadata <- expand.grid(
  p = 50,
  theta_1 = c(0),
  theta_2 = c(0, 0.1),
  theta_3 = c(0),
  nclust = c(1, 4),
  n = 50000
) %>% mutate(scenario = paste0("S", 1:4)) %>% dplyr::select(-nclust)

# -----------------------------
# Load and process all results
# -----------------------------
all_files <- list.files(results_dir, pattern = paste0("^", scenario_prefix, ".*\\.rds$"), full.names = TRUE)

read_and_expand <- function(path) {
  df <- readRDS(path)
  df$source_file <- path
  return(df)
}

all_results <- map_dfr(all_files, read_and_expand) %>%
  mutate(scenario = str_extract(source_file, "S\\d+(_null)?")) %>%
  left_join(scenario_metadata, by = "scenario") %>%
  select(-source_file) %>%
  mutate(
    ebic = if_else(
      str_detect(method, "^res_fused_lasso_"),
      str_replace(method, "^res_fused_lasso_", "") %>%
        str_replace("^0(\\d)$", "0.\\1") %>%   # turn 05 into 0.5, 01 into 0.1, etc.
        as.numeric(),
      NA_real_
    ),
    method = if_else(
      str_detect(method, "^res_fused_lasso_"),
      "res_fused_lasso",
      method
    )
  )

# -----------------------------
# Summary statistics
# -----------------------------
summary_stats_theta_2 <- all_results %>%
  group_by(scenario, method, ebic) %>%
  summarise(
    theta_2 = mean(theta_2),
    mean = mean(Estimate, na.rm = TRUE),
    bias = mean(Estimate - theta_2, na.rm = TRUE),
    sd = sd(Estimate, na.rm = TRUE),
    coverage = mean(theta_2 >= LCL & theta_2 <= UCL, na.rm = TRUE),
    ci_width = mean(UCL - LCL, na.rm = TRUE),
    nclust = mean(nclust, na.rm = TRUE),
    rand = mean(rand, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

filter(summary_stats_theta_2, scenario == "S1") %>%
    print(width = Inf)

filter(summary_stats_theta_2, scenario == "S2") %>%
  print(width = Inf)

filter(summary_stats_theta_2, scenario == "S3") %>%
  print(width = Inf)

filter(summary_stats_theta_2, scenario == "S4") %>%
  print(width = Inf)

# -----------------------------
# Plots
# -----------------------------

## 1. Boxplot of theta_hat by method
ggplot(all_results, aes(x = method, y = theta_hat, fill = method)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_hline(aes(yintercept = theta_true), color = "red", linetype = "dashed") +
  facet_wrap(~ scenario, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = expression(hat(theta)~"distribution across methods"),
       y = expression(hat(theta)),
       x = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # ggsave(file.path(plot_dir, "boxplot_theta_hat_all_methods.png"), width = 12, height = 6)

## 2. Bias and RMSE
summary_long <- summary_stats %>%
  pivot_longer(cols = c(bias, rmse), names_to = "metric", values_to = "value")

ggplot(filter(summary_long, scenario %in% c("S3", "S4")), aes(x = method, y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(metric ~ scenario, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "Bias and RMSE per method",
       x = "Method", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # +
  # ggsave(file.path(plot_dir, "bias_rmse_all_methods.png"), width = 12, height = 6)

# -----------------------------
# Done
# -----------------------------
message("✅ Summary and plots saved in: ", results_dir, " and ", plot_dir)
