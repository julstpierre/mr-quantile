#-------------------------------------------------------------------------------------------------------------------
# R script to run our STAN implementation of MR-Horse: Bayesian Lasso for Mendelian Randomization as described in:
# Grant, AJ and Burgess, S (2024). A Bayesian approach to Mendelian randomization 
# using summary statistics in the univariable and multivariable settings with correlated pleiotropy. 
# American Journal of Human Genetics. 111(1):165-180. doi: https://doi.org/10.1016/j.ajhg.2023.12.002
#-------------------------------------------------------------------------------------------------------------------

#-----------------------
# Load libraries
#-----------------------
library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)

#-----------------------
# Installing CmdStan
#-----------------------
check_cmdstan_toolchain()
cmdstan_path()
cmdstan_version()

#-----------------------
# Compiling the model
# ----------------------
file <- file.path(getwd(), "MR-Horse.stan")
mod <- cmdstan_model(file)
mod$print()
mod$exe_file()

#-----------------------
# Simulate data
#-----------------------
# Simulation parameters
set.seed(123)
args <- data.frame(frac = 0.25, beta = 0.2)
beta <- args$beta
n <- 10000
p <- 30
phi <- rep(0, p) 
alpha <- rnorm(p, 0, 0.15)
valid.IVs <- sample(1:p, (1-args$frac)*p)
invalid.IVs <- setdiff(1:p, valid.IVs)
phi[valid.IVs] <- 0
alpha[valid.IVs] <- 0

# Simulate genetic data and confounder
maf <- runif(p, 0.1, 0.5)
G <- sapply(1:p, function(j) rbinom(2 * n, 2, maf[j]))
U <- G %*% phi + rnorm(2 * n)

# Simulate uniform genetic effects such that the total proportion of variance explained 
# in the exposure by direct effects of the genetic variants is on average 10%
sigma2_gamma <- 1/9 * (var(U) + 1) * n / sum(G^2)
gamma <- runif(p, 0, sqrt(12 * sigma2_gamma))

# Simulate exposure and outcome
X <- G %*% gamma + U + rnorm(2 * n)
Y <- G %*% alpha + beta * X + U + rnorm(2 * n)

# Summary statistics
sampX = seq(1, n)
sampY = seq((n+1), (2*n))

# Summary stats for X
beta_hat_x <- unname(sapply(1:p, function(i) coef(lm(X[sampX] ~ G[sampX, i]))[-1]))
sigma2_x <- unname(sapply(1:p, function(i) coefficients(summary(lm(X[sampX] ~ G[sampX, i])))[2, "Std. Error"]^2))

# Summary stats for Y
beta_hat_y <- unname(sapply(1:p, function(i) coef(lm(Y[sampY] ~ G[sampY, i]))[-1]))
sigma2_y <- unname(sapply(1:p, function(i) coefficients(summary(lm(Y[sampY] ~ G[sampY, i])))[2, "Std. Error"]^2))

# Fit model
data_list <- list(p= p, 
                  Gamma_hat_y = beta_hat_y, 
                  gamma_hat_x = beta_hat_x, 
                  sigma_y = sqrt(sigma2_y),
                  sigma_x = sqrt(sigma2_x),
                  scale_global = 1,
                  slab_scale = 0.3,
                  slab_df = 4
                  )
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.9,
  init = rep(list(list(beta = median(beta_hat_y / beta_hat_x), alpha = rep(0, p))), 4),
  refresh = 500 # print update every 500 iters,
)

#----------------------------
# Posterior summary statistics
#----------------------------
# Summary of the posterior distribution
fit$summary(
  variables = c("beta", "alpha", "tau", "c"),
  "mean", "sd", "median",
  pep = ~ 2 * pmin(mean(. < 0), mean(. > 0)),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
) %>% print(n=Inf)

# Sort alpha's by posterior error prob. and select invalid instruments with a FDR < 0.05
fit$summary(variables = "alpha", "mean", pep = ~ 2 * pmin(mean(. < 0), mean(. > 0))) %>%
  arrange(pep) %>%
  mutate(qvalue = cummean(pep)) %>%
  filter(qvalue < .05) %>%
  mutate(index = as.integer(gsub("alpha\\[(\\d+)\\]", "\\1", variable))) %>%
  arrange(index) %>%
  select(-pep, -qvalue, -index)

# Compare with true non-zero alpha's
data.frame(variable = paste0("alpha[", which(alpha != 0), "]"), mean = round(alpha[which(alpha != 0)], 4)) %>%
  tibble()

# Plot of number of invalid instruments as a function of the FDR
library(ggplot2)
fit$summary(variables = "alpha", pep = ~ 2 * pmin(mean(. < 0), mean(. > 0))) %>%
  arrange(pep) %>%
  mutate(qvalue = cummean(pep),  n = row_number()) %>%
  ggplot(aes(x = qvalue, y = n)) +
  geom_line() +
  labs(x = "FDR", y = "Number of invalid instruments") +
  theme_minimal()

#----------------------------
# Posterior draws
#----------------------------
# Extracting draws
draws_df <- fit$draws(format = "df")

# Trace plots
mcmc_trace(fit$draws("alpha"))
mcmc_trace(fit$draws("beta"))

# Histograms
mcmc_hist(fit$draws("beta"))
mcmc_hist(fit$draws("alpha"))

#----------------------------
# Sampler diagnostics
#----------------------------
fit$sampler_diagnostics()
fit$diagnostic_summary()
