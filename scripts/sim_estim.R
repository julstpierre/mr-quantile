#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(MendelianRandomization)
library(quantreg)
library(glmnet)
# library(mrclust)
library(genlasso)
library(fossil)
library(mclust)
library(posterior)

# Load custom estimator functions
source("utils.R")
source("mrclust.R")
source("mrlaplace.R")

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript sim_estim.R <scenario_name> <start_chunk> <end_chunk>")
}

scenario_name <- args[1]          # e.g., "D_S1"
start_idx     <- as.numeric(args[2])  # e.g., 1
end_idx       <- as.numeric(args[3])  # e.g., 200

# -----------------------------
# Load simulation data
# -----------------------------
for (f in paste0("D_S", 1:8, ".RData")) {
  load(paste0("../data/simulations/", f))
  rm(f)
}
D <- get(scenario_name)                # Extract the dataset corresponding to scenario_name

# Safety check for indexing
end_idx <- min(end_idx, length(D))
chunk <- D[start_idx:end_idx]

# -----------------------------
# Run all MR methods
# -----------------------------
res_ivw         <- lapply(chunk, theta_ivw) %>% do.call(rbind.data.frame, .)
res_median      <- lapply(chunk, theta_median) %>% do.call(rbind.data.frame, .)
res_egger      <- lapply(chunk, theta_egger) %>% do.call(rbind.data.frame, .)
res_mrclust     <- lapply(chunk, fit_mrclust) %>% do.call(rbind.data.frame, .)
res_mrlaplace <- lapply(chunk, fit_mrlaplace) %>% do.call(rbind.data.frame, .)
res_mrlaplace_null <- lapply(chunk, fit_mrlaplace, null_cluster = TRUE) %>% do.call(rbind.data.frame, .)
res_mrlaplace_bayes <- lapply(chunk, fit_mrlaplace_bayes) %>% do.call(rbind.data.frame, .)

# -----------------------------
# Combine results into a single data frame and add a column for the model
# -----------------------------
tag_method <- function(res_list, method) {
  if (length(res_list) <= 8) {
    bind_rows(res_list) %>%
      mutate(method = method)
  } else {
    # Bind all elements first
    df <- bind_rows(res_list)
    
    # Split into two parts: first 4 cols and last 4 cols
    df_main <- df[, -c(5:8), drop = FALSE] %>%
      mutate(method = method)
    
    df_median <- df[, -c(1:4), drop = FALSE] %>%
      mutate(method = paste0(method, "_median"))
    
    bind_rows(df_main, df_median)
  }
}

res_all <- bind_rows(
  tag_method(res_ivw, "ivw")
  ,tag_method(res_median, "median")
  ,tag_method(res_egger, "egger")
  ,tag_method(res_mrclust, "mrclust")
  ,tag_method(res_mrlaplace, "mrlaplace")
  ,tag_method(res_mrlaplace_null, "mrlaplace_null")
  ,tag_method(res_mrlaplace_bayes, "mrlaplace_bayes")
  )

# -----------------------------
# Save results
# -----------------------------
out_dir <- "../results/rds"
out_file <- paste0(scenario_name, "_chunk_", start_idx, "_to_", end_idx, ".rds")
out_path <- file.path(out_dir, out_file)
saveRDS(res_all, out_path)

message("✅ Done: ", scenario_name, " | chunk ", start_idx, "–", end_idx, " | saved to ", out_path)