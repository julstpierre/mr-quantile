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
source("../../R/mr_quantile.R")

# Define tagging function to add method name to results
tag_method <- function(res_list, method) {
  bind_rows(res_list) %>% 
    mutate(method = method)
}

# -----------------------------
# Parse command-line arguments
# -----------------------------
scenario <- "weak_invalid_IVs"

if (scenario == "main"){
  scenario_ds <- paste0("D_S", 1:42)
} else if (scenario == "weak_invalid_IVs"){
  scenario_ds <- paste0("D_S", 1:56)
}

for (scenario_name in scenario_ds){

  # -----------------------------
  # Load simulation data
  # -----------------------------
  load(paste0("../data/", scenario, "/",scenario_name, ".RData"))
  chunk <- get(scenario_name)
  rm(list = scenario_name)
  
  # -----------------------------
  # Run all MR methods
  # -----------------------------
  res_all <- tag_method(lapply(chunk, fit_mr_quantile_vm), "mr_quantile_vm")
  
  # -----------------------------
  # Save results
  # -----------------------------
  out_dir <- paste0("../results/rds/", scenario)
  out_file <- paste0(scenario_name, "_mr_quantile_vm.rds")
  out_path <- file.path(out_dir, out_file)
  saveRDS(res_all, out_path)
  
  message("✅ Done: ", scenario_name, " | saved to ", out_path)
  
}