library(parallel)
library(MASS)
library(dplyr)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 6
num_iters <- 1000
params <- expand.grid(
  p = c(30, 100),
  beta = c(0, 0.4),
  h2 = c(0.30),
  nclust = c(1, 4),
  n = 50000,
  theta_M = 0
)

#----------------------------------------------------------------------
# Simulation scenarios
#----------------------------------------------------------------------

set.seed(2025070101)
D_S1_null <- filter(params, beta == 0, nclust == 1, p == 30) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025070102)
D_S1 <- filter(params, beta == 0.4, nclust == 1, p == 30) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025070103)
D_S2_null <- filter(params, beta == 0, nclust == 1, p == 100) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025070104)
D_S2 <- filter(params, beta == 0.4, nclust == 1, p == 100) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025070105)
D_S3 <- filter(params, beta == 0.4, nclust == 4, p == 30) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025070106)
D_S4 <- filter(params, beta == 0.4, nclust == 4, p == 100) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

# Save simulated data
save(D_S1, D_S1_null, D_S2, D_S2_null, D_S3, D_S4, file = '../data/simulations/sim_dat.R')
