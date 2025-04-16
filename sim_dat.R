library(parallel)
library(MASS)
library(dplyr)
source("functions.R")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
num_jobs <- 6
num_iters <- 1000
params <- expand.grid(
  frac = c(0.1, 0.3, 0.5, 0.7),
  scenario = 1:3,
  theta1 = c(0, 0.2)
)

#----------------------------------------------------------------------
# Theta = 0
#----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")

# Scenario 1
# 10% invalid IVs
set.seed(2024112501)
D_S1_10_null <- filter(params, scenario == 1, frac == 0.1, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112502)
D_S1_30_null <- filter(params, scenario == 1, frac == 0.3, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112503)
D_S1_50_null <- filter(params, scenario == 1, frac == 0.5, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112504)
D_S1_70_null <- filter(params, scenario == 1, frac == 0.7, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Scenario 2
# 10% invalid IVs
set.seed(2024112505)
D_S2_10_null <- filter(params, scenario == 2, frac == 0.1, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112506)
D_S2_30_null <- filter(params, scenario == 2, frac == 0.3, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112507)
D_S2_50_null <- filter(params, scenario == 2, frac == 0.5, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112508)
D_S2_70_null <- filter(params, scenario == 2, frac == 0.7, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Scenario 3
# 10% invalid IVs
set.seed(2024112509)
D_S3_10_null <- filter(params, scenario == 3, frac == 0.1, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112510)
D_S3_30_null <- filter(params, scenario == 3, frac == 0.3, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112511)
D_S3_50_null <- filter(params, scenario == 3, frac == 0.5, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112512)
D_S3_70_null <- filter(params, scenario == 3, frac == 0.7, theta1 == 0) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

save(D_S1_10_null, D_S1_30_null, D_S1_50_null, D_S1_70_null, D_S2_10_null,
     D_S2_30_null, D_S2_50_null, D_S2_70_null, D_S3_10_null, D_S3_30_null,
     D_S3_50_null, D_S3_70_null, file = 'Data/sim_dat_null.R')

#----------------------------------------------------------------------
# Theta = 0.2
#----------------------------------------------------------------------
# Scenario 1
# 10% invalid IVs
set.seed(2024112513)
D_S1_10 <- filter(params, scenario == 1, frac == 0.1, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112514)
D_S1_30 <- filter(params, scenario == 1, frac == 0.3, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112515)
D_S1_50 <- filter(params, scenario == 1, frac == 0.5, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112516)
D_S1_70 <- filter(params, scenario == 1, frac == 0.7, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Scenario 2
# 10% invalid IVs
set.seed(2024112517)
D_S2_10 <- filter(params, scenario == 2, frac == 0.1, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112518)
D_S2_30 <- filter(params, scenario == 2, frac == 0.3, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112519)
D_S2_50 <- filter(params, scenario == 2, frac == 0.5, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112520)
D_S2_70 <- filter(params, scenario == 2, frac == 0.7, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Scenario 3
# 10% invalid IVs
set.seed(2024112521)
D_S3_10 <- filter(params, scenario == 3, frac == 0.1, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 30% invalid IVs
set.seed(2024112522)
D_S3_30 <- filter(params, scenario == 3, frac == 0.3, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 50% invalid IVs
set.seed(2024112523)
D_S3_50 <- filter(params, scenario == 3, frac == 0.5, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# 70% invalid IVs
set.seed(2024112524)
D_S3_70 <- filter(params, scenario == 3, frac == 0.7, theta1 == 0.2) %>%
                list() %>%
                replicate(num_iters, .) %>%
                mclapply(mc.cores = num_jobs, sim_dat, mc.set.seed = TRUE)

# Save simulated data
save(D_S1_10, D_S1_30, D_S1_50, D_S1_70, D_S2_10, D_S2_30, D_S2_50, D_S2_70,
     D_S3_10, D_S3_30, D_S3_50, D_S3_70, file = 'Data/sim_dat.R')
