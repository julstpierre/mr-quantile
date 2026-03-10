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
  df <- readRDS(path)
  df$source_file <- path
  
  df <- df %>%
    mutate(
      across(c(Pvalue, Estimate, StdError, LCL, UCL),
             ~ suppressWarnings(as.numeric(.)))
    )
  
  return(df)
}

all_results <- map_dfr(all_files, read_and_expand) %>%
  mutate(scenario = str_extract(source_file, "S\\d+(_null)?")) %>%
  left_join(scenario_metadata, by = "scenario") %>%
  select(-source_file) %>%
  filter(!method %in% c("mr_quantreg_ma"))

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
  
  filter(summary_stats, p == 30, n == 50000, q == 0, InSIDE == TRUE) %>%
    print(width = Inf, n = Inf)
  
  filter(summary_stats, p == 30, n == 50000, q == 0.2) %>%
    print(width = Inf, n = Inf)
  
  filter(summary_stats, p == 30, n == 50000, q == 0.4) %>%
    print(width = Inf, n = Inf)
  
  filter(summary_stats, p == 30, n == 50000, q == 0.6, InSIDE == FALSE) %>%
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

## 2. Plot Power
plot_power <- function(
    data,
    theta_val = c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2),
    h2_u_val = NULL,
    n_val,
    p_val,
    methods = NULL,
    facet_row,
    facet_col,
    facet_labels = list(),
    title = ""
) {
  
  facet_formula <- as.formula(paste(facet_row, "~", facet_col))
  
  title_expr <- if (is.null(h2_u_val)) {
    bquote("Power, n = " * .(n_val) * ", " * p == .(p_val)  * .(title))
  } else {
    bquote("Power, n = " * .(n_val) * ", " * h[u]^2 == .(h2_u_val) * .(title))
  }
  
  ggplot(
    dplyr::filter(
      data,
      theta %in% theta_val,
      n %in% n_val,
      p %in% p_val,
      if (!is.null(methods)) method %in% methods else TRUE,
      if (!is.null(h2_u_val)) h2_u %in% h2_u_val else TRUE
    ) %>%
        mutate(method = factor(method, levels = methods)),
    aes(x = theta, y = power, fill = method)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    geom_text(
      aes(label = round(power, 2)),
      position = position_dodge(width = 0.9),
      vjust = -0.3,
      size = 2.5
    ) +
    coord_cartesian(ylim = c(0, 1.1)) +
    facet_grid(
      facet_formula,
      labeller = labeller(!!!facet_labels)
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = title_expr,
      y = "Percentage significant",
      x = "Method"
    ) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

if (sim_name == "weak_invalid_IVs"){
    plot_power(
      data = summary_stats,
      h2_u_val = 0,
      n = 50000,
      p = 50,
      facet_row = "h2_y",
      facet_col = "theta",
      facet_labels = list(
        h2_y = function(x) paste0("h²ᵧ = ", x),
        theta = function(x) paste0("θ = ", x)
      ),
      title = ", InSIDE Hold"
    )
    
    plot_power(
      data = summary_stats,
      h2_u_val = 0.1,
      n = 50000,
      p = 50,
      facet_row = "h2_y",
      facet_col = "theta",
      facet_labels = list(
        h2_y = function(x) paste0("h²ᵧ = ", x),
        theta = function(x) paste0("θ = ", x)
      ),
      title = ", InSIDE Violated"
    )
} else if (sim_name == "main"){
    plot_power(
      data = filter(summary_stats, InSIDE == TRUE),
      n = 50000,
      p = 30,
      facet_row = "q",
      facet_col = "theta",
      facet_labels = list(
        theta = function(x) paste0("θ = ", x),
        q = function(x) paste0("Invalid IVs = ", 100 * as.numeric(x), "%")
      ),
      title = ", InSIDE Hold"
    )
  
    plot_power(
      data = filter(summary_stats, InSIDE == FALSE),
      n = 50000,
      p = 30,
      facet_row = "q",
      facet_col = "theta",
      facet_labels = list(
        theta = function(x) paste0("θ = ", x),
        q = function(x) paste0("Invalid IVs = ", 100 * as.numeric(x), "%")
      ),
      title = ", InSIDE Violated"
   )
}

## 3. Plot Coverage
plot_coverage <- function(
    data,
    theta_val = c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2),
    h2_u_val = NULL,
    n_val,
    p_val = 50,
    methods = NULL,
    facet_row,
    facet_col,
    facet_labels = list(),
    title= ""
) {
  
  facet_formula <- as.formula(paste(facet_row, "~", facet_col))
  
  title_expr <- if (is.null(h2_u_val)) {
    bquote("Coverage, n = " * .(n_val) * ", " * p == .(p_val)  * .(title))
  } else {
    bquote("Coverage, n = " * .(n_val) * ", " * h[u]^2 == .(h2_u_val) * .(title))
  }
  
  ggplot(
    dplyr::filter(
      data,
      theta %in% theta_val,
      n %in% n_val,
      p %in% p_val,
      if (!is.null(methods)) method %in% methods else TRUE,
      if (!is.null(h2_u_val)) h2_u %in% h2_u_val else TRUE
    ),
    aes(x = theta, y = coverage, fill = method)
  ) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    geom_text(
      aes(label = round(coverage, 2)),
      position = position_dodge(width = 0.9),
      vjust = -0.3,
      size = 2.5
    ) +
    coord_cartesian(ylim = c(0.85, 1.05)) +
    # geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_grid(
      facet_formula,
      labeller = labeller(!!!facet_labels)
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = title_expr,
      y = "Coverage",
      x = "Method"
    ) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

if (sim_name == "weak_invalid_IVs"){
  plot_coverage(
    data = summary_stats,
    h2_u_val = 0,
    n = 50000,
    facet_row = "h2_y",
    facet_col = "theta",
    facet_labels = list(
      h2_y = function(x) paste0("h²ᵧ = ", x),
      theta = function(x) paste0("θ = ", x)
    )
  )
  
  plot_coverage(
    data = summary_stats,
    h2_u_val = 0.1,
    n = 50000,
    facet_row = "h2_y",
    facet_col = "theta",
    facet_labels = list(
      h2_y = function(x) paste0("h²ᵧ = ", x),
      theta = function(x) paste0("θ = ", x)
    )
  )
} else if (sim_name == "main"){
  
  plot_coverage(
    data = filter(summary_stats, InSIDE == TRUE),
    n = 50000,
    p = 30,
    facet_row = "q",
    facet_col = "theta",
    facet_labels = list(
      theta = function(x) paste0("θ = ", x),
      q = function(x) paste0("Invalid IVs = ", 100 * as.numeric(x), "%")
    ),
    title = ", InSIDE Hold"
  )
  
  plot_coverage(
    data = filter(summary_stats, InSIDE == FALSE),
    n = 50000,
    p = 30,
    facet_row = "q",
    facet_col = "theta",
    facet_labels = list(
      theta = function(x) paste0("θ = ", x),
      q = function(x) paste0("Invalid IVs = ", 100 * as.numeric(x), "%")
    ),
    title = ", InSIDE Violated"
  )
}

# Plot bias
library(ggpubr)

# Function to create a single plot
create_plot <- function(filter_expr, title = "") {
  # Filter data for this combination
  plot_data <- all_results %>%
    filter({{ filter_expr }}) %>%
    filter(method %in% method_levels) %>%
    mutate(bias = Estimate - theta) %>%
    add_row(method = "",) %>%
    mutate(method = factor(method, levels = c("", method_levels)))
  
  # Calculate statistics
  label_df <- plot_data %>%
    group_by(method) %>%
    summarise(
      mean_bias = mean(bias, na.rm = TRUE),
      sd = sd(Estimate, na.rm = TRUE),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      power = mean(LCL > 0 | UCL < 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create axis labels
  bias_axis_labels <- setNames(
    paste0(round(label_df$mean_bias, 3), "\n", 
           round(label_df$sd, 3), "\n", 
           round(label_df$rmse, 3),"\n", 
           round(label_df$power, 3)),
    label_df$method
  )
  bias_axis_labels[1] <- "Bias\nSD\nRMSE\nPower "
  
  # Create plot
  p <- plot_data %>%
    ggplot(aes(x = method, y = Estimate, fill = method, color = method)) +
    geom_jitter(size = 0.2, alpha = 0.5, width = 0.1) +
    geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA, color = "black") +
    scale_x_discrete(labels = bias_axis_labels) +
    geom_hline(aes(yintercept = theta), linetype = "dashed", color = "black") +
    coord_cartesian(ylim = c(plot_data$theta[1] - 0.1, plot_data$theta[1] + 0.1)) +
    scale_fill_manual(values = method_colors, labels = method_labels, drop = FALSE) +
    scale_color_manual(values = method_colors, labels = method_labels, drop = FALSE) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      y = "",
      x = ""
    ) +
    theme(
      axis.text.x = element_text(hjust = 0, size = 6),
      legend.position = "none",
      legend.title = element_blank(),
      plot.title = element_text(
        hjust = 0.5,   # 0 = left, 0.5 = center, 1 = right
        size  = 11,    # adjust to taste
        face  = "plain"
      )
    )
  
  return(p)
}

# Main simulations
p1 <- create_plot(filter_expr = p == 30 & q == 0.6 & InSIDE == FALSE & theta ==  0.1, title = expression(theta[0] == 0.1 * ", InSIDE Violated"))
p2 <- create_plot(filter_expr = p == 30 & q == 0.6 & InSIDE == TRUE & theta ==  0.1, title = expression(theta[0] == 0.1 * ", InSIDE Hold"))
p3 <- create_plot(filter_expr = p == 30 & q == 0.6 & InSIDE == FALSE & theta ==  -0.1, title = expression(theta[0] == -0.1 * ", InSIDE Violated"))
p4 <- create_plot(filter_expr = p == 30 & q == 0.6 & InSIDE == TRUE & theta ==  -0.1, title = expression(theta[0] == -0.1 * ", InSIDE Hold"))

# Combine plots
combined_plot <- ggpubr::ggarrange(
  p1, p2, p3, p4,
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

# Add overall title if needed
annotate_figure(combined_plot,
                bottom = text_grob("Method", size = 11, vjust = -10, hjust = -0.4),
                left = text_grob("Estimate", size = 11, rot = 90, hjust = -0.5, vjust = 2)
)

# Save plot
ggsave(filename = paste0(plot_dir, "/boxplots_main.pdf"),
               width = 10,
               height = 7
      )

# Weak invalid IVs
p1 <- create_plot(filter_expr = theta == -0.2 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == -0.2 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
                  )

p2 <- create_plot(filter_expr = theta == 0.2 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == 0.2 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
)

p3 <- create_plot(filter_expr = theta == -0.1 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == -0.1 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
)

p4 <- create_plot(filter_expr = theta == 0.1 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == 0.1 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
)

p5 <- create_plot(filter_expr = theta == -0.05 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == -0.05 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
)

p6 <- create_plot(filter_expr = theta == 0.05 & h2_u == 0.1 & h2_y == 0.1, 
                  title = expression(theta[0] == 0.05 ~ "," ~ h[y]^2 == 0.1 ~ "," ~ h[u]^2 == 0.1)
)
            
# Combine plots
combined_plot <- ggpubr::ggarrange(
  p1, p2, p3, p4, p5, p6,
  ncol = 2, nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# Add overall title if needed
annotate_figure(combined_plot,
                bottom = text_grob("Method", size = 11, vjust = -10, hjust = -0.4),
                left = text_grob("Estimate", size = 11, rot = 90, hjust = -0.5, vjust = 2)
)

# Save plot
ggsave(filename = paste0(plot_dir, "/boxplots_weak_invalid_IVs.pdf"),
       width = 10,
       height = 7
)

# -----------------------------
# Done
# -----------------------------
message("✅ Summary and plots saved in: ", results_dir, " and ", plot_dir)
