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
num_jobs <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "6"))
num_iters <- 5000
params <- expand.grid(
  p = 50,
  theta_1 = c(0.2),
  theta_2 = c(0, 0.8),
  theta_3 = c(-1, 1),
  nclust = c(1, 4),
  n = 50000,
  h2 = 0.10
)

#----------------------------------------------------------------------
# Simulation scenarios
#----------------------------------------------------------------------
scenario_list <- list(
  filter(params, theta_2 == 0,   theta_3 == -1, nclust == 1), # S1
  filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 1), # S2
  filter(params, theta_2 == 0,   theta_3 == 1,  nclust == 1), # S3
  filter(params, theta_2 == 0.8, theta_3 == 1,  nclust == 1), # S4
  filter(params, theta_2 == 0,   theta_3 == -1, nclust == 4), # S5
  filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 4), # S6
  filter(params, theta_2 == 0,   theta_3 == 1,  nclust == 4), # S7
  filter(params, theta_2 == 0.8, theta_3 == 1,  nclust == 4)  # S8
)

df_param <- scenario_list[[scenario_id]] %>% list()

set.seed(2025081300 + scenario_id)

res <- replicate(num_iters, df_param) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

# Save file
D_name <- paste0("D_S", scenario_id)
assign(paste0("D_S", scenario_id), res)
save(list = D_name, file = sprintf("../data/simulations/D_S%d.RData", scenario_id))
