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
library(ggpubr)

getwd()

# -----------------------------
# Parameters
# -----------------------------
sim_name <- c("main", "weak_invalid_IVs")[2]
results_dir <- paste0("../results/rds/", sim_name)
plot_dir <- paste0("../results/plots/", sim_name)
scenario_prefix <- "D_S"

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# Scenario metadata
# -----------------------------
if (sim_name == "weak_invalid_IVs"){
  scenario_metadata <- expand.grid(
    p = 50, 
    n = 50000, 
    theta = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2), 
    h2_x = 0.5, 
    h2_y = c(0.1, 0.2, 0.4, 0.6), 
    h2_u = c(0, 0.1), 
    q = 0.6
  ) %>% mutate(scenario = paste0("S", row_number()))
  
  group_vars <- c("theta", "p", "n", "h2_y", "h2_u", "method")
} else if (sim_name == "main"){
  scenario_metadata <- expand.grid(
    p = c(30, 100),
    theta = c(-0.1, 0, 0.1), 
    n = 50000,
    q = c(0, 0.2, 0.4, 0.6),
    InSIDE = c(TRUE, FALSE)
  ) %>%
    filter(!(InSIDE == FALSE & q == 0)) %>% 
    mutate(scenario = paste0("S", row_number()))
  
  group_vars <- c("InSIDE","theta", "p", "n", "q", "method")
}

# -----------------------------
# Load and process all results
# -----------------------------
all_files <- list.files(results_dir, pattern = paste0("^", scenario_prefix, ".*\\.rds$"), full.names = TRUE)

read_and_expand <- function(path) {
  readRDS(path) %>%
    filter(!method %in% c("mr_quantreg_ma", "mr_quantile_vm")) %>%
    mutate(
      across(any_of(c("Pvalue", "Estimate", "StdError", "LCL", "UCL")),
             ~ suppressWarnings(as.numeric(.))),
      scenario = str_extract(basename(path), "S\\d+(_null)?")
    )
}

all_results <- map(all_files, read_and_expand, .progress = TRUE) %>%
  list_rbind() %>%
  left_join(scenario_metadata, by = "scenario", relationship = "many-to-one")

# -----------------------------
# Summary statistics
# -----------------------------
summary_stats <- all_results %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    theta = mean(theta),
    mean = mean(Estimate, na.rm = TRUE),
    bias = mean(Estimate - theta, na.rm = TRUE),
    se = mean(StdError, na.rm = TRUE),
    sd = sd(Estimate, na.rm = TRUE),
    coverage = mean(theta >= LCL & theta <= UCL, na.rm = TRUE),
    power = if (all(is.na(Pvalue))) {
      mean(LCL > 0 | UCL < 0, na.rm = TRUE)
    } else {
      mean(Pvalue < 0.05, na.rm = TRUE)
    },
    ci_width = mean(UCL - LCL, na.rm = TRUE),
    iter = n(),
    .groups = "drop"
  )

if (sim_name == "weak_invalid_IVs") {
  
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
  
} else if (sim_name == "main") {
  
  filter(summary_stats, n == 50000, q == 0) %>%
    filter(method %in% c("mr_quantile_vm", "mr_quantreg")) %>% 
    print(width = Inf, n = Inf) 
  
  filter(summary_stats, n == 50000, q == 0.2, InSIDE == TRUE) %>%
    filter(method %in% c("mr_quantile_vm", "mr_quantreg")) %>% 
    print(width = Inf, n = Inf)
  
  filter(summary_stats, n == 50000, q == 0.2, InSIDE == FALSE) %>%
    filter(method %in% c("mr_quantile_vm", "mr_quantreg")) %>% 
    print(width = Inf, n = Inf)
  
  filter(summary_stats, p == 30, n == 50000, q == 0.4) %>%
    filter(method %in% c("mr_quantile_vm", "mr_quantreg")) %>% 
    print(width = Inf, n = Inf)
  
  filter(summary_stats, p == 30, n == 50000, q == 0.6) %>%
    filter(method %in% c("mr_quantile_vm", "mr_quantreg")) %>% 
    print(width = Inf, n = Inf)
}

# -----------------------------
# Plots
# -----------------------------

## 0. Create color palette
method_levels <- c(
  "mr_quantreg", "mr_cML", "mr_cML_DP", "mr_contmix", "mr_egger", "mr_ivw",
  "mr_ivw_oracle", "mr_lasso", "mr_median", "mr_mix",
  "mr_mode", "mr_presso", "mr_raps"
)

method_colors <- c(
  mr_cML        = "#1F78B4",
  mr_cML_DP     = "#6A3D9A",
  mr_contmix    = "#2F4F4F",
  mr_egger      = "#E31A1C",
  mr_ivw        = "#FF7F00",
  mr_ivw_oracle = "grey",
  mr_lasso      = "#A6CEE3",
  mr_median     = "#B2DF8A",
  mr_mix        = "#00BFC4",
  mr_mode       = "#FDBF6F",
  mr_presso     = "#33A02C",
  mr_quantreg   = "#CC79A7",
  mr_quantile_vm = "#F8766D",
  mr_raps       = "#B15928"
)

method_labels <- c(
  mr_cML        = "cML-MA-BIC",
  mr_cML_DP     = "cML-MA-BIC-DP",
  mr_contmix    = "MR-ContMix",
  mr_egger      = "MR-Egger",
  mr_ivw        = "MR-IVW",
  mr_ivw_oracle = "MR-IVW (Oracle)",
  mr_lasso      = "MR-Lasso",
  mr_median     = "MR-Weighted-Median",
  mr_mix        = "MR-Mix",
  mr_mode       = "MR-Weighted-Mode",
  mr_presso     = "MR-PRESSO",
  mr_quantreg   = "MR-Quantile",
  mr_quantile_vm   = "MR-Quantile-VM",
  mr_raps       = "MR-RAPS"
)

## 1. Plot Type-I error
plot_typeI_error <- function(
    data, 
    theta_val = 0, 
    n_val,
    facet_row,
    facet_col,
    facet_labels = list(),
    title = "",
    methods = method_levels
) {
  
  facet_formula <- as.formula(paste(facet_row, "~", facet_col))
  
  ggplot(filter(data, theta %in% theta_val, n %in% n_val, method %in% methods) %>%
           mutate(method = factor(method, levels = methods)), 
         aes(x = theta, y = power, fill = method)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(aes(yintercept = 0.05), linetype = "dashed", color = "black") +
    scale_fill_manual(
      values = method_colors, 
      labels = method_labels, 
      drop = FALSE
    ) +
    geom_text(
      aes(
        y = -Inf,                 # anchor to bottom of panel
        label = round(power, 3)
      ),
      position = position_dodge(width = 0.9),
      vjust = 0.5,               # push text further down
      angle = 45,                 # rotate text (change as you like)
      size = 2
    ) +
    coord_cartesian(clip = "off") +
    facet_grid(
      facet_formula,
      labeller = labeller(!!!facet_labels)
    ) +
    theme_minimal(base_size = 11) +
    labs(
      # title = bquote("Type-I error, " * theta[0] == .(theta_val) * .(title)),
      y = "% Type-I error",
      x = "Method"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 15, unit = "pt")),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(r = 10, unit = "pt", l = 5))
    )
}

if (sim_name == "weak_invalid_IVs"){
  plot_typeI_error(
    summary_stats,
    n=50000,
    facet_row = "h2_u",
    facet_col = "h2_y",
    facet_labels = list(
      h2_y = function(x) paste0("h²y = ", x),
      h2_u = function(x) paste0("h²u = ", x)
    )
  )
  
  # ggsave(filename = paste0(plot_dir, "/typeI_error.pdf"),
  #        width = 10,
  #        height = 7
  # )
  
} else if (sim_name == "main"){
  plot_typeI_error(
    filter(summary_stats, InSIDE == TRUE),
    n = 50000,
    facet_row = "p",
    facet_col = "q",
    facet_labels = list(
      p = function(x) paste0("p = ", x),
      q = function(x) paste0(100 * as.numeric(x), "% ", "Invalid IVs")
    )
  )
  
  # ggsave(filename = paste0(plot_dir, "/typeI_error_InSIDE_hold.pdf"),
  #        width = 10,
  #        height = 7
  # )
  
  plot_typeI_error(
    filter(summary_stats, InSIDE == FALSE),
    n = 50000,
    facet_row = "p",
    facet_col = "q",
    facet_labels = list(
      p = function(x) paste0("p = ", x),
      q = function(x) paste0(100 * as.numeric(x), "% ", "Invalid IVs")
    )
  )
  
  # ggsave(filename = paste0(plot_dir, "/typeI_error_InSIDE_violated.pdf"),
  #        width = 10, 
  #        height = 7
  # )
}

# Plot bias
library(patchwork)
library(tidyr)

metric_levels <- c("Bias", "SD", "RMSE", "Power", "Coverage")

# 3 decimals without the leading zero (".025"), and no "-0.000"
fmt <- function(v) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", round(v, 3) + 0))

create_plot <- function(filter_expr, title = "", text_size = 2.8) {
  plot_data <- all_results %>%
    filter({{ filter_expr }}, method %in% method_levels) %>%
    mutate(bias = Estimate - theta,
           method = factor(method, levels = method_levels))
  
  theta0 <- unique(plot_data$theta)
  stopifnot(length(theta0) == 1)
  
  stats_df <- plot_data %>%
    group_by(method) %>%
    summarise(
      Bias     = mean(bias, na.rm = TRUE),
      SD       = sd(Estimate, na.rm = TRUE),
      RMSE     = sqrt(mean(bias^2, na.rm = TRUE)),
      Power    = mean(LCL > 0 | UCL < 0, na.rm = TRUE),
      Coverage = mean(LCL < theta & UCL > theta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-method, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = rev(metric_levels)),
           label  = fmt(value))
  
  p_box <- ggplot(plot_data, aes(method, Estimate, fill = method, color = method)) +
    geom_jitter(size = 0.2, alpha = 0.5, width = 0.1) +
    geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA, color = "black") +
    geom_hline(yintercept = theta0, linetype = "dashed") +
    coord_cartesian(ylim = theta0 + c(-0.1, 0.1)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = method_colors, labels = method_labels,
                      drop = FALSE, name = NULL) +
    scale_color_manual(values = method_colors, labels = method_labels,
                       drop = FALSE, name = NULL) +
    guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
    labs(title = title, x = NULL, y = "Estimate") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 11))
  
  p_tab <- ggplot(stats_df, aes(method, metric, label = label)) +
    geom_text(size = text_size) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 9, colour = "black"))
  
  p_box / p_tab + plot_layout(heights = c(3, 1.2))
}

p1 <- create_plot(p == 30 & q == 0.6 & InSIDE == FALSE & theta ==  0.1,
                  expression(theta[0] == 0.1 * ", InSIDE Violated"))
p2 <- create_plot(p == 30 & q == 0.6 & InSIDE == TRUE  & theta ==  0.1,
                  expression(theta[0] == 0.1 * ", InSIDE Holds"))
p3 <- create_plot(p == 30 & q == 0.6 & InSIDE == FALSE & theta == -0.1,
                  expression(theta[0] == -0.1 * ", InSIDE Violated"))
p4 <- create_plot(p == 30 & q == 0.6 & InSIDE == TRUE  & theta == -0.1,
                  expression(theta[0] == -0.1 * ", InSIDE Holds"))

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "boxplots_main.pdf"), combined_plot,
       width = 10, height = 8)

# Weak invalid IVs
theta_values <- c(-0.2, 0.2, -0.1, 0.1, -0.05, 0.05)

plots_weak <- lapply(theta_values, function(th) {
  create_plot(
    abs(theta - th) < 1e-8 & h2_u == 0.1 & h2_y == 0.1,
    title = bquote(theta[0] == .(th) ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1),
    text_size = 2.8
  )
})

combined_plot <- wrap_plots(plots_weak, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11))

ggsave(file.path(plot_dir, "boxplots_weak_invalid_IVs.pdf"), combined_plot,
       width = 10, height = 11)

# -----------------------------
# Done
# -----------------------------
message("✅ Summary and plots saved in: ", results_dir, " and ", plot_dir)
