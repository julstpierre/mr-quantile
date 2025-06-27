
#' Fit BayesFuseMR: Bayesian fused shrinkage for Mendelian Randomization
#'
#' This function fits a Bayesian fused shrinkage model for Mendelian Randomization using
#' the cmdstanr interface. The input summary statistics are sorted by scaled ratio estimates,
#' and a regularized horseshoe prior is applied.
#'
#' @param beta_hat_y Numeric vector of estimated SNP effects on the outcome.
#' @param beta_hat_x Numeric vector of estimated SNP effects on the exposure.
#' @param sigma_y Numeric vector of standard errors for beta_hat_y.
#' @param sigma_x Numeric vector of standard errors for beta_hat_x.
#' @param snps Numeric vector of standard errors for beta_hat_x. Default is paste0("SNP", 1:length(beta_hat_y)).
#' @param seed Optional integer. Random seed for reproducibility of the MCMC sampling.
#' @param scale_local Numeric scalar. Scale parameter for the prior on the first \code{theta} element. Default is 1.
#' @param scale_global Numeric scalar. Global scale parameter for the half-Cauchy prior on \code{tau} (overall shrinkage). Default is 1.
#' @param slab_scale Numeric scalar. Slab scale parameter for the regularized horseshoe prior. Default is 0.1.
#' @param slab_df Numeric scalar. Degrees of freedom for the slab t-distribution in the regularized horseshoe prior. Default is 4.
#'
#' @return A \code{cmdstanr} fitted model object containing posterior samples.
#'
#' @export
fit_BayesFuseMR <- function(beta_hat_y, beta_hat_x, sigma_y, sigma_x,
                            snps = paste0("SNP", 1:length(beta_hat_y)),
                            seed = NULL,
                            scale_local = 1,
                            scale_global = 1,
                            slab_scale = 0.1,
                            slab_df = 4){
  #-----------------------
  # Compiling the model
  # ----------------------
  stan_file <- system.file("stan", "BayesFuseMR.stan", package = "BayesFuseMR")
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  #-------------------------------------------------------
  # Sort the data by increasing value of the scaled ratios
  # -----------------------------------------------------
  ratios <- beta_hat_y / beta_hat_x
  ratios.se <- sigma_y / abs(beta_hat_x)
  scaled.ratios = ratios/ratios.se
  ord <- order(scaled.ratios)
  
  #-------------------------------------------------------
  # Run MCMC
  #-------------------------------------------------------
  data_list <- list(p = length(ratios), 
                    pi_scaled = scaled.ratios[ord], 
                    sigma_pi = ratios.se[ord], 
                    scale_local = scale_local,     # scale parameter for theta_1
                    scale_global = scale_global,   # scale for the half cauchy prior for tau
                    slab_scale = slab_scale,       # slab scale for the regularized horseshoe
                    slab_df = slab_df              # slab degrees of freedom for the regularized horseshoe
  )
  fit <- mod$sample(
    data = data_list,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.85,
    refresh = 500 # print update every 500 iters
  )
  
  return(list(posterior = fit, snps = snps[ord]))
}
