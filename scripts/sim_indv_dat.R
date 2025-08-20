#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Need scenario number 1-8 as argument")
scenario_id <- as.integer(args[1])

library(parallel)
library(MASS)
library(dplyr)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
num_iters <- 1000
params <- expand.grid(
  p = 60,
  theta_1 = c(0),
  theta_2 = c(0, 0.1),
  theta_3 = c(0),
  q = c(0, 0.2, 0.4, 0.6),
  n = 50000
)

#----------------------------------------------------------------------
# Simulation scenarios
#----------------------------------------------------------------------
#for (scenario_id in 1:8){
  df_param <- params[scenario_id, ] %>% list()
  
  set.seed(2025081300 + scenario_id)
  
  res <- replicate(num_iters, df_param) %>%
    mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)
  
  # Save file
  D_name <- paste0("D_S", scenario_id)
  assign(paste0("D_S", scenario_id), res)
  
  save(list = D_name, file = sprintf("../data/simulations/D_S%d.RData", scenario_id))
#}
