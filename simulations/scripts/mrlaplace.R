
#' Fit BayesClustMR: Bayesian Mendelian Randomization
#'
#' This function fits a Bayesian fused shrinkage model for Mendelian Randomization using
#' the cmdstanr interface. The input summary statistics are ratio estimates, and a mixture
#' of Laplace distributions is applied.
#'
#' @param beta_hat_y Numeric vector of estimated SNP effects on the outcome.
#' @param beta_hat_x Numeric vector of estimated SNP effects on the exposure.
#' @param sigma_y Numeric vector of standard errors for beta_hat_y.
#' @param sigma_x Numeric vector of standard errors for beta_hat_x.
#' @param snps Numeric vector of standard errors for beta_hat_x. Default is paste0("SNP", 1:length(beta_hat_y)).
#' @param seed Optional integer. Random seed for reproducibility of the MCMC sampling.
#'
#' @return A \code{cmdstanr} fitted model object containing posterior samples.
#'
#' @export
mr_laplace <- function(beta_hat_y, beta_hat_x, sigma2_y,
                            snps = paste0("SNP", 1:length(beta_hat_y)),
                            min_clust_search = 10,
                            stop_waic_iter = 5,
                            seed = NULL){
  #-----------------------
  # Compiling the models
  # ----------------------
  # stan_file <- system.file("stan", "BayesClustMR.stan", package = "BayesFuseMR")
  stan_file <- "../inst/stan/BayesClustMR.stan"
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  #-------------------------------------------------------
  # Prepare data for Stan
  #-------------------------------------------------------
  ratios <- beta_hat_y / beta_hat_x
  ratios.se <- sqrt( sigma2_y / beta_hat_x^2 )
  data_list <- list(N = length(ratios),
                    theta = ratios,
                    b = ratios.se / sqrt(2)
  )
  
  #-------------------------------------------------------
  # Fit the model over a range of K values
  #-------------------------------------------------------
  K_list <- 1:length(beta_hat_y)
  
  for (K in K_list){
    
    # Update number of clusters
    data_list$K <- K
    
    # Initial cluster assignments
    cluster <- cluster::pam(ratios, k = K, diss = FALSE, metric = "manhattan", cluster.only=TRUE)
    mu <- numeric(K)
    for (k in 1:K) {
      xk <- ratios[cluster == k]
      mu[k] <- weighted_quantile(xk, rep(1, length(xk)), prob = 0.5)
    }
    
    # Initial values
    init_list <- list(list(mu = sort(mu), 
                           pi = table(match(cluster, order(mu))) / length(ratios),
                           lambda = ratios.se
                           )
                      ) %>% replicate(4, .)
    
    # Run MCMC
    fit <- mod$sample(
      data = data_list,
      init = init_list,
      seed = seed,
      chains = 4,
      parallel_chains = 4,
      adapt_delta = 0.9,
      refresh = 500 # print update every 500 iters
    )
    
    # Extract log-likelihood
    log_lik <- fit$draws("log_lik", format = "matrix")
    
    # Compute WAIC
    waic_vals <- loo::waic(log_lik)
    
    # Create table with WAIC results
    if (K==1){
      best_waic <- waic_vals
      best_fit <- fit
      best_K <- 1
      waic_list <- data.frame(K = K,
                              elpd_waic = waic_vals$estimates["elpd_waic", "Estimate"],
                              elpd_diff = 0,
                              se_diff = 0,
                              best_K = 1
                              )
    # Compare WAIC with best model so far
    } else {
      waic_diff <- loo::loo_compare(list("best_model" = best_waic, "new_model" = waic_vals))
      if (rownames(waic_diff)[1] == "new_model") {waic_diff["best_model", "elpd_diff"] <- -waic_diff["best_model", "elpd_diff"]}
      waic_list <- rbind(waic_list, c(K=K, waic_vals$estimates["elpd_waic", "Estimate"], waic_diff[2, c(1, 2)], best_K = best_K))
      
      if (waic_list[K, "elpd_diff"] - 1.96 * waic_list[K, "se_diff"] > 0){
        best_waic <- waic_vals
        best_fit <- fit
        best_K <- K
      }
    }
    
    # Early stopping if no improvement in WAIC
    if (K >= min_clust_search && best_K != K){
      last_waic <- waic_list[(nrow(waic_list)-stop_waic_iter+1):nrow(waic_list), "best_K"]
      if (all(last_waic == last_waic[1])){
        message(paste0("No improvement in WAIC over last ", stop_waic_iter, " iterations. Stopping search."))
        break
      }
    }
  }
  
  # Summary of the posterior distribution
  posterior_summary <- best_fit$summary(
                        variables = c("mu", "pi"), "mean", "median", "sd", "rhat",
                        pep = ~ 2 * pmin(mean(. < 0), mean(. > 0)),
                        extra_quantiles = ~posterior::quantile2(., probs = c(.025, .975))
  ) 
  
  # Assign cluster labels based on responsibilities
  r_summary <- summarise_draws(best_fit$draws(variables = "r"), mean)
  r_matrix <- matrix(r_summary$mean, nrow = length(ratios), ncol = best_K, byrow = FALSE)
  cluster <- apply(r_matrix, 1, which.max)
  
  return(list(
    results = data.frame(
        SNP = snps,
        cluster = cluster,
        pi_hat = round(as.numeric(table(cluster) / length(snps)), 3)[cluster], 
        probability = round(apply(r_matrix, 1, max), 3),
        cluster_median = round(posterior_summary[cluster, "mean", drop = TRUE], 3),
        cluster_se = round(posterior_summary[cluster, "sd", drop = TRUE], 3),
        cluster_lcl = round(posterior_summary[cluster, "q2.5", drop = TRUE], 3),
        cluster_ucl = round(posterior_summary[cluster, "q97.5", drop = TRUE], 3),
        theta = ratios,
        theta_se = ratios.se
    ),
    posterior_summary = posterior_summary,
    waic_list = waic_list
  ))
  
}
