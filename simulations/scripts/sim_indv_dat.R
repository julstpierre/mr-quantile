#!/usr/bin/env Rscript

#=========================
# Load libraries & utils
#=========================
rm(list = ls())
library(dplyr)
source("utils.R")

#-----------------------------------
# Parse command-line argument
#-----------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide a scenario ID as the first argument.")

scenario_id <- as.integer(args[1])
cat("Running scenario D_S:", scenario_id, "\n")

#-----------------------------------
# Define the simulation parameters
#-----------------------------------
params <- expand.grid(
  p = c(30, 100),
  theta = c(-0.1, 0, 0.1), 
  n = 50000,
  q = c(0, 0.2, 0.4, 0.6),
  InSIDE = c(TRUE, FALSE)
) %>%
  filter(!(InSIDE == FALSE & q == 0)) %>%
  mutate(num_iter = ifelse(theta == 0, 10000, 2000))

df_param <- params[scenario_id, ] %>% list()

set.seed(2025121900 + scenario_id)

#----------------------------------------------------------------------
# Simulate data
#----------------------------------------------------------------------
res <- lapply(1:df_param[[1]]$num_iter, function(i) sim_indv_dat(df_param[[1]]))

#----------------------------------------------------------------------
# Save file
#----------------------------------------------------------------------
D_name <- paste0("D_S", scenario_id)
assign(paste0("D_S", scenario_id), res)
save(list = D_name, file = sprintf("../data/simulations/main/D_S%d.RData", scenario_id))
