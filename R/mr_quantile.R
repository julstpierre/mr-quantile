#' Mendelian randomization via quantile regression
#'
#' Estimates a causal effect using a quantile-based Mendelian randomization
#' estimator. The method operates on ratio estimates derived from SNP-exposure
#' and SNP-outcome summary statistics and iteratively estimates the causal
#' parameter using a weighted quantile regression framework.
#'
#' @param beta_hat_y Numeric vector of SNP-outcome association estimates.
#' @param beta_hat_x Numeric vector of SNP-exposure association estimates.
#' @param sigma2_y Numeric vector of variances of the SNP-outcome estimates.
#' @param sigma2_x Numeric vector of variances of the SNP-exposure estimates.
#' @param tau Optional quantile parameter. If `NULL`, it is estimated.
#' @param lambda Optional scale parameter. If `NULL`, it is estimated.
#' @param weights Optional vector of observation weights.
#' @param tol Convergence tolerance for the log-likelihood.
#' @param max_iter Maximum number of iterations for the optimization.
#' @param se Method for standard error estimation. Options are `"boot"` for
#' bootstrap or `NULL` for no standard error estimation.
#' @param nboot Number of bootstrap replicates when `se = "boot"`.
#' @param verbose Logical; if `TRUE`, prints optimization progress.
#'
#' @return A list containing:
#' \describe{
#'   \item{coefficients}{Data frame with the estimate, standard error,
#'   and 95\% confidence interval.}
#'   \item{tau}{Estimated quantile parameter.}
#'   \item{lambda}{Estimated scale parameter.}
#'   \item{aic}{Akaike information criterion.}
#'   \item{bic}{Bayesian information criterion.}
#' }
#'
#' @details
#' The method first computes ratio estimates from the summary statistics.
#' It then iteratively updates the quantile parameter, causal estimate,
#' and scale parameter by maximizing the asymmetric Laplace likelihood.
#'
#' Standard errors can be obtained via parametric bootstrap by sampling
#' new SNP-exposure and SNP-outcome associations from their estimated
#' sampling distributions.
#'
#' @examples
#' data("dummy_data")
#'
#' res <- mr_quantile(dummy_data$beta_hat_y, 
#'                    dummy_data$beta_hat_x, 
#'                    dummy_data$sigma_y^2, 
#'                    dummy_data$sigma_x^2,
#'                    verbose = TRUE
#')
#' res$coefficients
#'
#' @importFrom stats sd rnorm qnorm
#' @export
mr_quantile <- function(beta_hat_y,
                        beta_hat_x,
                        sigma2_y,
                        sigma2_x,
                        tau = NULL,
                        lambda = NULL,
                        weights = rep(1, length(beta_hat_y)),
                        tol = 1e-5,
                        max_iter = 500,
                        se = "boot",
                        nboot = 1000,
                        verbose = FALSE) {
  
  n <- length(beta_hat_y)
  
  # Compute ratios
  ratios <- beta_hat_y / beta_hat_x
  ratios.se <- sqrt(sigma2_y / beta_hat_x^2 +
                      beta_hat_y^2 * sigma2_x / beta_hat_x^4)
  w_i <- 1 / ratios.se
  
  # Initialize parameters
  lambda_init <- lambda
  lambda <- ifelse(is.null(lambda_init), 1, lambda_init)
  lambda_i <- lambda * w_i
  
  tau_init <- tau
  tau <- ifelse(is.null(tau_init), 0.5, tau_init)
  
  # Optimization loop
  loglik_prev <- -Inf
  
  for (iter in 1:max_iter) {
    
    # Update theta
    theta <- weighted_quantile(ratios, weights * w_i, tau)
    
    # Update lambda
    if (is.null(lambda_init)) {
      lambda <- sum(weights) /
        sum(weights * w_i * rho_tau(ratios - theta, tau))
      lambda_i <- lambda * w_i
    }
    
    # Update tau
    if (is.null(tau_init)) {
      a <- sum(weights * lambda_i * (ratios - theta))
      tau <- 0.5 - a / ( 2 * ( 2 * sum(weights) + sqrt ( a^2 + 4 * sum(weights)^2 ) ) )
      tau <- min(max(tau, 0.1), 0.9)
    }
    
    # Log-likelihood
    loglik <- sum(weights) * (log(lambda) + log(tau * (1 - tau))) -
      lambda * sum(weights * w_i * rho_tau(ratios - theta, tau))
    
    # Convergence
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  if (iter == max_iter)
    warning("Maximum iterations reached without convergence.")
  
  # Bootstrap SE
  if (is.null(se)) {
    theta_se <- NA
  } else if (se == "boot") {
    
    boot_mat <- replicate(nboot, {
      beta_hat_y_boot <- rnorm(n, beta_hat_y, sqrt(sigma2_y))
      beta_hat_x_boot <- rnorm(n, beta_hat_x, sqrt(sigma2_x))
      
      fit_b <- mr_quantile(beta_hat_y_boot,
                           beta_hat_x_boot,
                           sigma2_y,
                           sigma2_x,
                           tau = tau_init,
                           lambda = lambda_init,
                           se = NULL)
      
      c(theta = fit_b$coefficients$Estimate,
        tau   = fit_b$tau,
        lambda = fit_b$lambda)
    })
    
    theta_se <- sd(boot_mat["theta", ])
  }
  
  # Print results if verbose
  if (verbose) {
    cat("Converged after", iter, "iterations\n")
    cat("Estimated tau:", tau, "\n")
    cat("Estimated lambda:", lambda, "\n")
  }
  
  # Return list
  list(
    coefficients = data.frame(
      Estimate = theta,
      StdError = theta_se,
      LCL = theta - qnorm(0.975) * theta_se,
      UCL = theta + qnorm(0.975) * theta_se
    ),
    tau = tau,
    lambda = lambda,
    aic = -2 * loglik + 2 * (2 + is.null(tau_init)),
    bic = -2 * loglik + log(n) * (2 + is.null(tau_init))
  )
}

mr_quantile_vm <- function(beta_hat_y,
                           beta_hat_x,
                           sigma2_y,
                           sigma2_x,
                           tau = NULL,
                           weights = rep(1, length(beta_hat_y)),
                           tol = 1e-5,
                           max_iter = 500,
                           se = "boot",
                           nboot = 1000,
                           verbose = FALSE) {
  
  n <- length(beta_hat_y)
  
  # Compute ratios
  ratios <- beta_hat_y / beta_hat_x
  ratios.se <- sqrt(sigma2_y / beta_hat_x^2 + beta_hat_y^2 * sigma2_x / beta_hat_x^4)
  w_i <- 1 / ratios.se
  
  # Initialize parameters
  tau_init <- tau
  tau <- ifelse(is.null(tau_init), 0.5, tau_init)
  c_tau <- function(t) sqrt(1 - 2*t + 2*t^2) / (t * (1 - t))
  
  loglik_prev <- -Inf
  
  for (iter in 1:max_iter) {
    
    # Update lambda
    lambda_i <- c_tau(tau) * w_i
    
    # Update theta
    theta <- weighted_quantile(ratios, weights * lambda_i, tau)
    
    # Update tau
    if (is.null(tau_init)) {
      tau_obj <- function(t) {
        li <- c_tau(t) * w_i
        sum(weights * (log(li) + log(t * (1 - t)))) -
          sum(weights * li * rho_tau(ratios - theta, t))
      }
      opt <- optimize(tau_obj, interval = c(1e-4, 1 - 1e-4), maximum = TRUE)
      tau <- opt$maximum
    }
    
    # Log-likelihood
    loglik <- sum(weights * (log(lambda_i) + log(tau * (1 - tau)))) -
      sum(weights * lambda_i * rho_tau(ratios - theta, tau))
    
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  if (iter == max_iter) warning("Maximum iterations reached without convergence.")
  
  # Bootstrap SE
  if (is.null(se)) {
    theta_se <- NA
  } else if (se == "boot") {
    
    boot_mat <- replicate(nboot, {
      beta_hat_y_boot <- rnorm(n, beta_hat_y, sqrt(sigma2_y))
      beta_hat_x_boot <- rnorm(n, beta_hat_x, sqrt(sigma2_x))
      
      fit_b <- mr_quantile_vm(beta_hat_y_boot,
                           beta_hat_x_boot,
                           sigma2_y,
                           sigma2_x,
                           tau = tau_init,
                           se = NULL)
      
      c(theta = fit_b$coefficients$Estimate,
        tau   = fit_b$tau,
        lambda = fit_b$lambda)
    })
    
    theta_se <- sd(boot_mat["theta", ])
  }
  
  # Print results if verbose
  if (verbose) {
    cat("Converged after", iter, "iterations\n")
    cat("Estimated tau:", tau, "\n")
    cat("Estimated lambda:", c_tau(tau), "\n")
  }
  
  # Return list
  list(
    coefficients = data.frame(
      Estimate = theta,
      StdError = theta_se,
      LCL = theta - qnorm(0.975) * theta_se,
      UCL = theta + qnorm(0.975) * theta_se
    ),
    tau = tau,
    lambda =  c_tau(tau),
    aic = -2 * loglik + 2 * (1 + is.null(tau_init)),
    bic = -2 * loglik + log(n) * (1 + is.null(tau_init))
  )
}