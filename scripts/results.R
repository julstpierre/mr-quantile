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
  n = 25000, 
  theta = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2), 
  h2_x = 0.5, 
  h2_y = c(0.1, 0.2, 0.4, 0.6), 
  h2_u = c(0, 0.1), 
  q = 0.6
) %>% mutate(scenario = paste0("S", row_number()))

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
summary_stats <- all_results %>%
  group_by(theta, h2_y, h2_u, method) %>%
  summarise(
    theta = mean(theta),
    mean = mean(Estimate, na.rm = TRUE),
    bias = mean(Estimate - theta, na.rm = TRUE),
    sd = sd(Estimate, na.rm = TRUE),
    coverage = mean(theta >= LCL & theta <= UCL, na.rm = TRUE),
    power = mean(LCL > 0 | UCL < 0, na.rm = TRUE),
    ci_width = mean(UCL - LCL, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

filter(summary_stats, h2_u == "0" & h2_y == "0.1") %>%
    print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0" & h2_y == "0.2") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0" & h2_y == "0.4") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0" & h2_y == "0.6") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0.1" & h2_y == "0.1") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0.1" & h2_y == "0.2") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0.1" & h2_y == "0.4") %>%
  print(width = Inf, n = Inf)

filter(summary_stats, h2_u == "0.1" & h2_y == "0.6") %>%
  print(width = Inf, n = Inf)

# -----------------------------
# Plots
# -----------------------------

## 1. Plot Type-I error and Power
plot_power <- function(data, h2_u_val = "0", n = 25000) {
  
  ggplot(filter(data, h2_u == h2_u_val), aes(x = theta, y = power, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(aes(yintercept = 0.05), linetype = "dashed", color = "red") +
    facet_grid(
      h2_y ~ theta,
      labeller = labeller(
        h2_y = function(x) paste0("h²ᵧ = ", x),
        theta = function(x) paste0("θ = ", x)
      )
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = bquote("Type-I error and Power, n = " * .(n) * ", " * h[u]^2 == .(h2_u_val)),
      y = "Percentage significant",
      x = "Method"
    ) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

plot_power(summary_stats, h2_u_val = "0", n = 25000)
plot_power(summary_stats, h2_u_val = "0.1", n = 25000)

ggplot(filter(summary_stats, h2_u == "0.1"), aes(x = theta, y = power, fill = method)) +
  geom_bar(stat="identity", position = "dodge") +
  geom_hline(aes(yintercept = 0.05), linetype = "dashed", color = "red") +
  facet_grid(h2_y ~ theta,
             labeller = labeller(
               h2_y = function(x) {
                 paste0("h²ᵧ = ", x)
               },   # or custom function / named vector
               theta = function(x) {
                 paste0("θ = ", x)
               }
             )) +
  theme_minimal(base_size = 13) +
  labs(title = bquote("Type-I error and Power, n = 25000, " * h[u]^2 == 0.1),
       y = "Percentage significant",
       x = "Method") +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
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
