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
scenario_metadata <- tibble(
  scenario = paste0("S", 1:9),
  theta_D_true = c(0, 0.4, 0, 0.4, 0, 0.4, 0.4, 0.4, 0.4),
  theta_T_true = c(-0.2, 0.2, 0, 0.4, 0.2, 0.6, 0.2, 0.4, 0.6)
)

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
  mutate(theta_T_lcl = theta_T_hat - 1.96 * theta_T_se, theta_T_ucl = theta_T_hat + 1.96 * theta_T_se) %>%
  mutate(theta_D_lcl = theta_D_hat - 1.96 * theta_D_se, theta_D_ucl = theta_D_hat + 1.96 * theta_D_se)

# -----------------------------
# Summary statistics
# -----------------------------
summary_stats_theta_T <- all_results %>%
  group_by(scenario, method) %>%
  summarise(
    theta_T = mean(theta_T_true),
    mean = mean(theta_T_hat, na.rm = TRUE),
    bias = mean(theta_T_hat - theta_T_true, na.rm = TRUE),
    sd = sd(theta_T_hat, na.rm = TRUE),
    coverage = mean(theta_T_true >= theta_T_lcl & theta_T_true <= theta_T_ucl, na.rm = TRUE),
    ci_width = mean(theta_T_ucl - theta_T_lcl, na.rm = TRUE),
    nclust = mean(nclust, na.rm = TRUE),
    rand = mean(rand, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

summary_stats_theta_D <- all_results %>%
  group_by(scenario, method) %>%
  summarise(
    theta_D = mean(theta_D_true),
    mean = mean(theta_D_hat, na.rm = TRUE),
    bias = mean(theta_D_hat - theta_D_true, na.rm = TRUE),
    sd = sd(theta_D_hat, na.rm = TRUE),
    coverage = mean(theta_D_true >= theta_D_lcl & theta_D_true <= theta_D_ucl, na.rm = TRUE),
    ci_width = mean(theta_D_ucl - theta_D_lcl, na.rm = TRUE),
    nclust = mean(nclust, na.rm = TRUE),
    rand = mean(rand, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

filter(summary_stats_theta_T, scenario == "S1") %>%
    print(width = Inf)

filter(summary_stats_theta_T, scenario == "S2") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S3") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S4") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S5") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S6") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S7") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S8") %>%
  print(width = Inf)

filter(summary_stats_theta_T, scenario == "S9") %>%
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
