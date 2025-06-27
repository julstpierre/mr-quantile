#' Summarize clusters from BayesFuseMR posterior samples
#'
#' This function computes posterior summary statistics and identifies clusters
#' based on significant breaks in the fused parameter \code{eta}. It assigns
#' cluster IDs to instruments (SNPs), computes weighted cluster means, standard errors,
#' confidence intervals, and returns a detailed data frame.
#'
#' Multiple hypothesis testing for successive differences is controlled at FDR < 0.20.
#'
#' @param fit A fitted \code{cmdstanr} model object returned by \code{fit_BayesFuseMR}.
#' @param FDR Numeric scalar. False discovery rate threshold for identifying significant breaks (default 0.3).
#'
#' @return A data frame with one row per instrument containing:
#' \itemize{
#'   \item \code{instrument}: Instrument (SNP) identifier.
#'   \item \code{cluster}: Cluster ID assigned based on significant breaks.
#'   \item \code{cluster_size}: Number of instruments in the cluster.
#'   \item \code{cluster_mean}: Weighted mean of \code{theta} in the cluster.
#'   \item \code{cluster_se}: Standard error of the cluster mean.
#'   \item \code{cluster_lcl}, \code{cluster_ucl}: 95\% confidence interval for the cluster mean.
#'   \item \code{theta_mean}, \code{theta_se}: Posterior mean and standard error for each instrument.
#' }
#'
#' @import dplyr
#' @importFrom posterior quantile2
#' @export
clust_BayesFuseMR <- function(fit, FDR = .3) {
  #-----------------------------
  # Posterior summary statistics
  #-----------------------------
  # Multiple hypothesis testing for successive difference at a FDR < 0.20
  breaks <-
    fit$posterior$summary(variables = "eta", "mean", 
                pep = ~ 2 * pmin(mean(. < 0), mean(. > 0)),
                extra_quantiles = ~posterior::quantile2(., probs = c(.25, .75))
                ) %>%
    arrange(pep) %>%
    mutate(qvalue = cummean(pep), ENFD = row_number() * qvalue) %>%
    filter(qvalue < FDR) %>%
    mutate(index = as.integer(gsub("eta\\[(\\d+)\\]", "\\1", variable))) %>%
    arrange(index) %>%
    pull(index)
  
  # Assign cluster IDs
  theta_df <- fit$posterior$summary(variables = "theta", "mean", "sd", pep = ~ 2 * pmin(mean(. < 0), mean(. > 0))) %>%
    rename(theta_mean = mean, theta_se = sd, theta_pep = pep) %>%
    mutate(snp = fit$snps) %>%
    mutate(cluster = cut(row_number(), breaks = c(0, breaks, max(row_number())), labels = FALSE, right = TRUE)) %>%
    select(snp, theta_mean, theta_se, theta_pep, cluster)
  
  # Compute cluster summary statistics
  theta_summary <-
    theta_df %>%
    mutate(weight = 1 / (theta_se^2)) %>%
    group_by(cluster) %>%
    summarise(
      cluster_mean = sum(theta_mean * weight) / sum(weight),
      cluster_se = sqrt(1 / sum(weight)),
      cluster_size = n(),
      .groups = "drop"
    ) %>%
    mutate(cluster_lcl = cluster_mean - 1.96 * cluster_se,
           cluster_ucl = cluster_mean + 1.96 * cluster_se) %>%
    left_join(theta_df, by = "cluster") %>%
    select(snp, cluster, cluster_size, cluster_mean, cluster_se, cluster_lcl, cluster_ucl, theta_mean, theta_se, theta_pep)
  
  return(theta_summary)
}
