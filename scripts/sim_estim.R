# Load libraries
library(dplyr)
library(quantreg)
library(glmnet)
library(MendelianRandomization)
library(MASS)
library(parallel)
library(mrclust)
library(genlasso)
library(BayesFuseMR)
library(fossil)
library(mclust)
source("utils.R")

# # Load data
load("../data/simulations/sim_dat.R")

#--------------------------------------------------------------------------------
# Scenario 1
#--------------------------------------------------------------------------------
res_S1_null_ivw <- mclapply(D_S1_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_null_egger <- mclapply(D_S1_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_null_median <- mclapply(D_S1_null, theta_median) %>% do.call(rbind, .)
res_S1_null_mrclust <- mclapply(D_S1_null, theta_mrclust) %>% do.call(rbind.data.frame, .)
res_S1_null_fused_lasso <- mclapply(D_S1_null, theta_fused_lasso) %>% do.call(rbind, .)
res_S1_null_BayesFuseMR <- mclapply(D_S1_null, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)

res_S1_ivw <- mclapply(D_S1, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_egger <- mclapply(D_S1, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_median <- mclapply(D_S1, theta_median) %>% do.call(rbind, .)
res_S1_mrclust <- mclapply(D_S1, theta_mrclust) %>% do.call(rbind.data.frame, .)
res_S1_fused_lasso <- mclapply(D_S1, theta_fused_lasso) %>% do.call(rbind, .)
res_S1_BayesFuseMR <- mclapply(D_S1, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)

#--------------------------------------------------------------------------------
# Scenario 2
#--------------------------------------------------------------------------------
res_S2_null_ivw <- mclapply(D_S2_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_null_egger <- mclapply(D_S2_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_null_median <- mclapply(D_S2_null, theta_median) %>% do.call(rbind, .)
res_S2_null_mrclust <- mclapply(D_S2_null, theta_mrclust) %>% do.call(rbind.data.frame, .)
res_S2_null_fused_lasso <- mclapply(D_S2_null, theta_fused_lasso) %>% do.call(rbind, .)
res_S2_null_BayesFuseMR <- mclapply(D_S2_null, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)

res_S2_ivw <- mclapply(D_S2, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_egger <- mclapply(D_S2, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_median <- mclapply(D_S2, theta_median) %>% do.call(rbind, .)
res_S2_mrclust <- mclapply(D_S2, theta_mrclust) %>% do.call(rbind.data.frame, .)
res_S2_fused_lasso <- mclapply(D_S2, theta_fused_lasso) %>% do.call(rbind, .)
res_S2_BayesFuseMR <- mclapply(D_S2, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)

#--------------------------------------------------------------------------------
# Scenario 3
#--------------------------------------------------------------------------------
res_S3_ivw <- mclapply(D_S3, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_egger <- mclapply(D_S3, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_median <- mclapply(D_S3, theta_median) %>% do.call(rbind, .)
res_S3_mrclust <- mclapply(D_S3, theta_mrclust) %>% do.call(rbind, .)
res_S3_fused_lasso <- mclapply(D_S3, theta_fused_lasso) %>% do.call(rbind, .)
res_S3_BayesFuseMR <- mclapply(D_S3, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)

#--------------------------------------------------------------------------------
# Scenario 4
#--------------------------------------------------------------------------------
res_S4_ivw <- mclapply(D_S4, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S4_egger <- mclapply(D_S4, theta_egger) %>% do.call(rbind.data.frame, .)
res_S4_median <- mclapply(D_S4, theta_median) %>% do.call(rbind, .)
res_S4_mrclust <- mclapply(D_S4, theta_mrclust) %>% do.call(rbind, .)
res_S4_fused_lasso <- mclapply(D_S4, theta_fused_lasso) %>% do.call(rbind, .)
res_S4_BayesFuseMR <- mclapply(D_S4, theta_BayesFuseMR, FDR = 0.05) %>% do.call(rbind, .)
