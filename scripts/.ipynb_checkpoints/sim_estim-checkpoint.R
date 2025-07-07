#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(MendelianRandomization)
library(quantreg)
library(glmnet)
library(mrclust)
library(genlasso)
library(BayesFuseMR)
library(fossil)
library(mclust)

# Load custom estimator functions (you should define these in utils.R)
source("utils.R")

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript sim_estim.R <scenario_name> <start_chunk> <end_chunk>")
}

scenario_name <- args[1]          # e.g., "D_S1_null"
start_idx     <- as.numeric(args[2])  # e.g., 1
end_idx       <- as.numeric(args[3])  # e.g., 200

# -----------------------------
# Load simulation data
# -----------------------------
load("../data/simulations/sim_dat.R")  # Loads D_S1_null, D_S1, ..., D_S4
D <- get(scenario_name)                # Extract the dataset corresponding to scenario_name

# Safety check for indexing
end_idx <- min(end_idx, length(D))
chunk <- D[start_idx:end_idx]

# -----------------------------
# Run all MR methods
# -----------------------------
res_ivw         <- lapply(chunk, theta_ivw) %>% do.call(rbind.data.frame, .)
res_egger       <- lapply(chunk, theta_egger) %>% do.call(rbind.data.frame, .)
res_median      <- lapply(chunk, theta_median) %>% do.call(rbind, .) %>% data.frame()
res_mrclust     <- lapply(chunk, theta_mrclust) %>% do.call(rbind.data.frame, .)
res_fused_lasso <- lapply(chunk, theta_fused_lasso) %>% do.call(rbind.data.frame, .)
res_bayesfusemr   <- lapply(chunk, theta_BayesFuseMR, FDR = c(0.05, 0.1, 0.2, 0.3, 0.4)) %>% do.call(rbind.data.frame, .)

names(res_median) <- names(res_ivw)

# -----------------------------
# Combine results into a single data frame and add a column for the model
# -----------------------------
tag_method <- function(res_list, method) {
  bind_rows(res_list) %>%
    mutate(method = method)
}

res_all <- bind_rows(
  tag_method(res_ivw, "ivw"),
  tag_method(res_egger, "egger"),
  tag_method(res_median, "median"),
  tag_method(res_mrclust, "mrclust"),
  tag_method(res_fused_lasso, "fused_lasso"),
  tag_method(res_bayesfusemr, "bayesfusemr")
)

# -----------------------------
# Save results
# -----------------------------
out_dir <- "../results"
out_file <- paste0(scenario_name, "_chunk_", start_idx, "_to_", end_idx, ".rds")
out_path <- file.path(out_dir, out_file)
saveRDS(res_all, out_path)

message("✅ Done: ", scenario_name, " | chunk ", start_idx, "–", end_idx, " | saved to ", out_path)