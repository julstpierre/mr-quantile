#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)
library(quantreg)
library(fossil)
library(mclust)
library(MRcML)

# Load custom estimator functions
source("utils.R")
# source("mrclust.R")

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript sim_estim.R <scenario> <scenario_name> <start_chunk> <end_chunk>")
}

scenario <- args[1]               # e.g., "main"
scenario_name <- args[2]          # e.g., "D_S1"
start_idx     <- as.numeric(args[3])  # e.g., 1
end_idx       <- as.numeric(args[4])  # e.g., 200

# -----------------------------
# Load simulation data
# -----------------------------
for (f in list.files(paste0("../data/simulations/", scenario, "/"), pattern = "^D_S[0-9]+\\.RData$")) {
  load(paste0("../data/simulations/", scenario, "/", f))
  rm(f)
}
D <- get(scenario_name)                # Extract the dataset corresponding to scenario_name

# Safety check for indexing
end_idx <- min(end_idx, length(D))
chunk <- D[start_idx:end_idx]

# -----------------------------
# Run all MR methods
# -----------------------------
res_mr_ivw        <- lapply(chunk, fit_mr_ivw) %>% do.call(rbind.data.frame, .)
res_mr_ivw_oracle <- lapply(chunk, fit_mr_ivw, snps = 31:50) %>% do.call(rbind.data.frame, .)
res_mr_median     <- lapply(chunk, fit_mr_median) %>% do.call(rbind.data.frame, .)
res_mr_egger      <- lapply(chunk, fit_mr_egger) %>% do.call(rbind.data.frame, .)
res_mr_mode       <- lapply(chunk, fit_mr_mode) %>% do.call(rbind.data.frame, .)
res_mr_raps       <- lapply(chunk, fit_mr_raps) %>% do.call(rbind.data.frame, .)
res_mr_lasso      <- lapply(chunk, fit_mr_lasso) %>% do.call(rbind.data.frame, .)
res_mr_contmix    <- lapply(chunk, fit_mr_contmix) %>% do.call(rbind.data.frame, .)
res_mr_mix        <- lapply(chunk, fit_mr_mix) %>% do.call(rbind.data.frame, .)
res_mr_quantreg_wlbs   <- lapply(chunk, fit_mr_quantreg, se = "wlbs") %>% do.call(rbind.data.frame, .)
res_mr_quantreg_boot   <- lapply(chunk, fit_mr_quantreg, se = "boot") %>% do.call(rbind.data.frame, .)
res_mr_quantreg_norm   <- lapply(chunk, fit_mr_quantreg, se = NULL) %>% do.call(rbind.data.frame, .)
res_mr_cML        <- lapply(chunk, fit_mr_cML) %>% do.call(rbind.data.frame, .)
res_mr_cML_DP     <- lapply(chunk, fit_mr_cML_DP) %>% do.call(rbind.data.frame, .)
res_mr_presso <- lapply(chunk, fit_mr_presso) %>% do.call(rbind.data.frame, .)

# ----------------------------
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
   tag_method(res_mr_ivw, "mr_ivw")
  ,tag_method(res_mr_ivw_oracle, "mr_ivw_oracle")
  ,tag_method(res_mr_median, "mr_median")
  ,tag_method(res_mr_egger, "mr_egger")
  ,tag_method(res_mr_mode, "mr_mode")
  ,tag_method(res_mr_raps, "mr_raps")
  ,tag_method(res_mr_lasso, "mr_lasso")
  ,tag_method(res_mr_contmix, "mr_contmix")
  ,tag_method(res_mr_mix, "mr_mix")
  ,tag_method(res_mr_quantreg_wlbs, "mr_quantreg_wlbs")
  ,tag_method(res_mr_quantreg_boot, "mr_quantreg_boot")
  ,tag_method(res_mr_quantreg_norm, "mr_quantreg_norm")
  ,tag_method(res_mr_cML, "mr_cML")
  ,tag_method(res_mr_cML_DP, "mr_cML_DP")
  ,tag_method(res_mr_presso, "mr_presso")
  )

# -----------------------------
# Save results
# -----------------------------
out_dir <- paste0("../results/rds/", scenario)
out_file <- paste0(scenario_name, "_chunk_", start_idx, "_to_", end_idx, ".rds")
out_path <- file.path(out_dir, out_file)
saveRDS(res_all, out_path)

message("✅ Done: ", scenario_name, " | chunk ", start_idx, "–", end_idx, " | saved to ", out_path)