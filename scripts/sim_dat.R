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
  p = c(50),
  theta_1 = c(0.2),
  theta_2 = c(0, 0.8),
  theta_3 = c(-1, 0, 1),
  h2 = c(0.30),
  nclust = c(1, 4),
  n = 50000
)

#----------------------------------------------------------------------
# Simulation scenarios
#----------------------------------------------------------------------

set.seed(2025070101)
D_S1 <- filter(params, theta_2 == 0, theta_3 == -1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070102)
D_S2 <- filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070103)
D_S3 <- filter(params, theta_2 == 0, theta_3 == 0, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070104)
D_S4 <- filter(params, theta_2 == 0.8, theta_3 == 0, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070105)
D_S5 <- filter(params, theta_2 == 0, theta_3 == 1, nclust == 1) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070106)
D_S6 <- filter(params, theta_2 == 0.8, theta_3 == 1, nclust == 1) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070107)
D_S7 <- filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070108)
D_S8 <- filter(params, theta_2 == 0.8, theta_3 == 0, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

set.seed(2025070109)
D_S9 <- filter(params, theta_2 == 0.8, theta_3 == 1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Save simulated data
save(D_S1, D_S2, D_S3, D_S4, D_S5, D_S6, D_S7, D_S8, D_S9, 
     file = '../data/simulations/sim_dat.R')
