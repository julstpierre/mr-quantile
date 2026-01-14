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
methods <- args[5]              # e.g., "{fast, slow}"

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
tag_method <- function(res_list, method) {
    bind_rows(res_list) %>% 
      mutate(method = method)
}

if (methods == "fast") {

  res_all <- bind_rows(
    tag_method(lapply(chunk, fit_mr_ivw), "mr_ivw"),
    tag_method(
      lapply(chunk, function(sim) {
        fit_mr_ivw(sim, snps = setdiff(1:sim$p, sim$invalid.IVs))
      }),
      "mr_ivw_oracle"
    ),
    tag_method(lapply(chunk, fit_mr_median), "mr_median"),
    tag_method(lapply(chunk, fit_mr_egger), "mr_egger"),
    tag_method(lapply(chunk, fit_mr_mode), "mr_mode"),
    tag_method(lapply(chunk, fit_mr_raps), "mr_raps"),
    tag_method(lapply(chunk, fit_mr_lasso), "mr_lasso"),
    tag_method(lapply(chunk, fit_mr_contmix), "mr_contmix"),
    tag_method(lapply(chunk, fit_mr_mix), "mr_mix"),
    tag_method(lapply(chunk, fit_mr_quantreg, se = "wlbs"), "mr_quantreg_wlbs"),
    tag_method(lapply(chunk, fit_mr_quantreg, se = "boot"), "mr_quantreg_boot"),
    tag_method(lapply(chunk, fit_mr_cML), "mr_cML")
  )
}

if (methods == "slow") {

  res_all <- bind_rows(
    tag_method(lapply(chunk, fit_mr_presso), "mr_presso"),
    tag_method(lapply(chunk, fit_mr_cML_DP), "mr_cML_DP")
  )
}

# -----------------------------
# Save results
# -----------------------------
out_dir <- paste0("../results/rds/", scenario)
out_file <- paste0(scenario_name, "_chunk_", start_idx, "_to_", end_idx, "_", methods, ".rds")
out_path <- file.path(out_dir, out_file)
saveRDS(res_all, out_path)

message("✅ Done: ", scenario_name, " | chunk ", start_idx, "–", end_idx, " | saved to ", out_path)
