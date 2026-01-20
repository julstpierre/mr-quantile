rm(list = ls())
library(parallel)
library(dplyr)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 4
params <- expand.grid(
  p = c(50, 100),
  theta = c(-0.1, -0.05, 0, 0.05, 0.1), 
  n = c(50000, 100000),
  q = c(0, 0.2, 0.4, 0.6)
) %>%
  mutate(num_iter = ifelse(theta == 0, 10000, 2000))

#----------------------------------------------------------------------
# Simulate data
#----------------------------------------------------------------------

for (scenario_id in 1:nrow(params)){
  df_param <- params[scenario_id, ] %>% list()
  
  set.seed(2025121900 + scenario_id)
  
  res <- replicate(df_param[[1]]$num_iter, df_param) %>%
    mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)
  
  # Save file
  D_name <- paste0("D_S", scenario_id)
  assign(paste0("D_S", scenario_id), res)
  save(list = D_name, file = sprintf("../data/simulations/main/D_S%d.RData", scenario_id))
}
