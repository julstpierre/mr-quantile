rm(list = ls())
library(parallel)
library(dplyr)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 4
num_iters <- 400
params <- expand.grid(
  p = 100000,
  theta = c(0, 0.2),
  n = 50000,
  h2_x = 0.3,
  h2_y = 0.3,
  m = 100,
  q = c(0, 0.2, 0.4, 0.6),
  pi = 0.2
)

#----------------------------------------------------------------------
# Simulate data
#----------------------------------------------------------------------

for (scenario_id in 1:nrow(params)){
  df_param <- params[scenario_id, ] %>% list()
  
  set.seed(2025082100 + scenario_id)
  
  res <- replicate(num_iters, df_param) %>%
    mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)
  
  # Save file
  D_name <- paste0("D_S", scenario_id)
  assign(paste0("D_S", scenario_id), res)
  save(list = D_name, file = sprintf("../data/simulations/D_S%d.RData", scenario_id))
}
