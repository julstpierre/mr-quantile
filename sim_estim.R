# Load libraries
library(dplyr)
library(quantreg)
library(glmnet)
library(MendelianRandomization)
library(MASS)
library(parallel)
library(hqreg)
source("functions.R")

# # Load data
load("Data/sim_dat_null.R")
load("Data/sim_dat.R")

#--------------------------------------------------------------------------------
# Scenario 1
#--------------------------------------------------------------------------------
# 10% invalid IVs
res_S1_10_null_ivw <- mclapply(D_S1_10_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_10_null_egger <- mclapply(D_S1_10_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_10_null_median <- mclapply(D_S1_10_null, theta_median) %>% do.call(rbind, .)
res_S1_10_null_lasso <- mclapply(D_S1_10_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S1_10_ivw <- mclapply(D_S1_10, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_10_egger <- mclapply(D_S1_10, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_10_median <- mclapply(D_S1_10, theta_median) %>% do.call(rbind, .)
res_S1_10_lasso <- mclapply(D_S1_10, theta_lasso) %>% do.call(rbind.data.frame, .)

# 30% invalid IVs
res_S1_30_null_ivw <- mclapply(D_S1_30_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_30_null_egger <- mclapply(D_S1_30_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_30_null_median <- mclapply(D_S1_30_null, theta_median) %>% do.call(rbind, .)
res_S1_30_null_lasso <- mclapply(D_S1_30_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S1_30_ivw <- mclapply(D_S1_30, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_30_egger <- mclapply(D_S1_30, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_30_median <- mclapply(D_S1_30, theta_median) %>% do.call(rbind, .)
res_S1_30_lasso <- mclapply(D_S1_30, theta_lasso) %>% do.call(rbind.data.frame, .)

# 50% invalid IVs
res_S1_50_null_ivw <- mclapply(D_S1_50_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_50_null_egger <- mclapply(D_S1_50_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_50_null_median <- mclapply(D_S1_50_null, theta_median) %>% do.call(rbind, .)
res_S1_50_null_lasso <- mclapply(D_S1_50_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S1_50_ivw <- mclapply(D_S1_50, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_50_egger <- mclapply(D_S1_50, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_50_median <- mclapply(D_S1_50, theta_median) %>% do.call(rbind, .)
res_S1_50_lasso <- mclapply(D_S1_50, theta_lasso) %>% do.call(rbind.data.frame, .)

# 70% invalid IVs
res_S1_70_null_ivw <- mclapply(D_S1_70_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_70_null_egger <- mclapply(D_S1_70_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_70_null_median <- mclapply(D_S1_70_null, theta_median) %>% do.call(rbind, .)
res_S1_70_null_lasso <- mclapply(D_S1_70_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S1_70_ivw <- mclapply(D_S1_70, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S1_70_egger <- mclapply(D_S1_70, theta_egger) %>% do.call(rbind.data.frame, .)
res_S1_70_median <- mclapply(D_S1_70, theta_median) %>% do.call(rbind, .)
res_S1_70_lasso <- mclapply(D_S1_70, theta_lasso) %>% do.call(rbind.data.frame, .)

#--------------------------------------------------------------------------------
# Scenario 2
#--------------------------------------------------------------------------------
# 10% invalid IVs
res_S2_10_null_ivw <- mclapply(D_S2_10_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_10_null_egger <- mclapply(D_S2_10_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_10_null_median <- mclapply(D_S2_10_null, theta_median) %>% do.call(rbind, .)
res_S2_10_null_lasso <- mclapply(D_S2_10_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S2_10_ivw <- mclapply(D_S2_10, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_10_egger <- mclapply(D_S2_10, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_10_median <- mclapply(D_S2_10, theta_median) %>% do.call(rbind, .)
res_S2_10_lasso <- mclapply(D_S2_10, theta_lasso) %>% do.call(rbind.data.frame, .)

# 30% invalid IVs
res_S2_30_null_ivw <- mclapply(D_S2_30_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_30_null_egger <- mclapply(D_S2_30_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_30_null_median <- mclapply(D_S2_30_null, theta_median) %>% do.call(rbind, .)
res_S2_30_null_lasso <- mclapply(D_S2_30_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S2_30_ivw <- mclapply(D_S2_30, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_30_egger <- mclapply(D_S2_30, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_30_median <- mclapply(D_S2_30, theta_median) %>% do.call(rbind, .)
res_S2_30_lasso <- mclapply(D_S2_30, theta_lasso) %>% do.call(rbind.data.frame, .)

# 50% invalid IVs
res_S2_50_null_ivw <- mclapply(D_S2_50_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_50_null_egger <- mclapply(D_S2_50_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_50_null_median <- mclapply(D_S2_50_null, theta_median) %>% do.call(rbind, .)
res_S2_50_null_lasso <- mclapply(D_S2_50_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S2_50_ivw <- mclapply(D_S2_50, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_50_egger <- mclapply(D_S2_50, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_50_median <- mclapply(D_S2_50, theta_median) %>% do.call(rbind, .)
res_S2_50_lasso <- mclapply(D_S2_50, theta_lasso) %>% do.call(rbind.data.frame, .)

# 70% invalid IVs
res_S2_70_null_ivw <- mclapply(D_S2_70_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_70_null_egger <- mclapply(D_S2_70_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_70_null_median <- mclapply(D_S2_70_null, theta_median) %>% do.call(rbind, .)
res_S2_70_null_lasso <- mclapply(D_S2_70_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S2_70_ivw <- mclapply(D_S2_70, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S2_70_egger <- mclapply(D_S2_70, theta_egger) %>% do.call(rbind.data.frame, .)
res_S2_70_median <- mclapply(D_S2_70, theta_median) %>% do.call(rbind, .)
res_S2_70_lasso <- mclapply(D_S2_70, theta_lasso) %>% do.call(rbind.data.frame, .)

#--------------------------------------------------------------------------------
# Scenario 3
#--------------------------------------------------------------------------------
# 10% invalid IVs
res_S3_10_null_ivw <- mclapply(D_S3_10_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_10_null_egger <- mclapply(D_S3_10_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_10_null_median <- mclapply(D_S3_10_null, theta_median) %>% do.call(rbind, .)
res_S3_10_null_lasso <- mclapply(D_S3_10_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S3_10_ivw <- mclapply(D_S3_10, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_10_egger <- mclapply(D_S3_10, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_10_median <- mclapply(D_S3_10, theta_median) %>% do.call(rbind, .)
res_S3_10_lasso <- mclapply(D_S3_10, theta_lasso) %>% do.call(rbind.data.frame, .)

# 30% invalid IVs
res_S3_30_null_ivw <- mclapply(D_S3_30_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_30_null_egger <- mclapply(D_S3_30_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_30_null_median <- mclapply(D_S3_30_null, theta_median) %>% do.call(rbind, .)
res_S3_30_null_lasso <- mclapply(D_S3_30_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S3_30_ivw <- mclapply(D_S3_30, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_30_egger <- mclapply(D_S3_30, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_30_median <- mclapply(D_S3_30, theta_median) %>% do.call(rbind, .)
res_S3_30_lasso <- mclapply(D_S3_30, theta_lasso) %>% do.call(rbind.data.frame, .)

# 50% invalid IVs
res_S3_50_null_ivw <- mclapply(D_S3_50_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_50_null_egger <- mclapply(D_S3_50_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_50_null_median <- mclapply(D_S3_50_null, theta_median) %>% do.call(rbind, .)
res_S3_50_null_lasso <- mclapply(D_S3_50_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S3_50_ivw <- mclapply(D_S3_50, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_50_egger <- mclapply(D_S3_50, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_50_median <- mclapply(D_S3_50, theta_median) %>% do.call(rbind, .)
res_S3_50_lasso <- mclapply(D_S3_50, theta_lasso) %>% do.call(rbind.data.frame, .)

# 70% invalid IVs
res_S3_70_null_ivw <- mclapply(D_S3_70_null, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_70_null_egger <- mclapply(D_S3_70_null, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_70_null_median <- mclapply(D_S3_70_null, theta_median) %>% do.call(rbind, .)
res_S3_70_null_lasso <- mclapply(D_S3_70_null, theta_lasso) %>% do.call(rbind.data.frame, .)

res_S3_70_ivw <- mclapply(D_S3_70, theta_ivw) %>% do.call(rbind.data.frame, .)
res_S3_70_egger <- mclapply(D_S3_70, theta_egger) %>% do.call(rbind.data.frame, .)
res_S3_70_median <- mclapply(D_S3_70, theta_median) %>% do.call(rbind, .)
res_S3_70_lasso <- mclapply(D_S3_70, theta_lasso) %>% do.call(rbind.data.frame, .)

#--------------------------------------------------------------------------------  
# Variable selection
#--------------------------------------------------------------------------------
# # Save results for selection of IVs
# lasso.IVs <- lasso.path.IVs[[lam_pos]]
# 
# # (a) The proportion of instruments selected as invalid which were truly invalid;
# precision <- length(intersect(setdiff(1:p, lasso.IVs), invalid.IVs)) / length(setdiff(1:p, lasso.IVs))
#   
# # (b) The proportion of instruments selected as valid which were truly valid;
# npv <- length(intersect(lasso.IVs, valid.IVs)) / length(lasso.IVs)
# 
# # (c) The proportion of truly invalid instruments which were selected as invalid;
# tpr <- length(intersect(setdiff(1:p, lasso.IVs), invalid.IVs)) / length(invalid.IVs)
# 
# # (d) The proportion of truly valid instruments which were estimated as valid.
# specificity <- length(intersect(lasso.IVs, valid.IVs)) / length(valid.IVs)
