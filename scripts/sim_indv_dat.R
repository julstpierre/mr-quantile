rm(list = ls())
library(parallel)
library(dplyr)
library(MASS)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 4
num_iters <- 200
params <- expand.grid(
  p = 90,
  theta_1 = c(0, 0.1),
  n = 10000,
  m = 10,
  k = 4,
  q = c(0, 0.2, 0.4, 0.6),
  rho = 0.5
)

#----------------------------------------------------------------------
# Simulate data
#----------------------------------------------------------------------

for (scenario_id in 1:nrow(params)){
  df_param <- params[scenario_id, ] %>% list()
  
  set.seed(2025082800 + scenario_id)
  
  res <- replicate(num_iters, df_param) %>%
    mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)
  
  # Save file
  D_name <- paste0("D_S", scenario_id)
  assign(paste0("D_S", scenario_id), res)
  save(list = D_name, file = sprintf("../data/simulations/D_S%d.RData", scenario_id))
}
