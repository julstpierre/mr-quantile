library(parallel)
library(MASS)
library(dplyr)
source("utils.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 6
num_iters <- 5000
params <- expand.grid(
  p = 50,
  theta_1 = c(0.2),
  theta_2 = c(0, 0.8),
  theta_3 = c(-1, 1),
  nclust = c(1, 4),
  n = 10000,
  h2 = 0.10
)

#----------------------------------------------------------------------
# Simulation scenarios
#----------------------------------------------------------------------

set.seed(2025081301)
D_S1 <- filter(params, theta_2 == 0, theta_3 == -1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081302)
D_S2 <- filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081303)
D_S3 <- filter(params, theta_2 == 0, theta_3 == 1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081304)
D_S4 <- filter(params, theta_2 == 0.8, theta_3 == 1, nclust == 1) %>%
          list() %>%
          replicate(num_iters, .) %>%
          mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081305)
D_S5 <- filter(params, theta_2 == 0, theta_3 == -1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081306)
D_S6 <- filter(params, theta_2 == 0.8, theta_3 == -1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081307)
D_S7 <- filter(params, theta_2 == 0, theta_3 == 1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

set.seed(2025081308)
D_S8 <- filter(params, theta_2 == 0.8, theta_3 == 1, nclust == 4) %>%
  list() %>%
  replicate(num_iters, .) %>%
  mclapply(mc.cores = num_jobs, sim_indv_dat, mc.set.seed = TRUE)

# Save simulated data
for (d in paste0("D_S", 1:8)) {
  save(list = d,
       file = paste0("../data/simulations/", d, ".RData")
       )
}
