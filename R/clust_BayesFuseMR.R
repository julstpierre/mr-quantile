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
#' @param FDR Numeric scalar. False discovery rate threshold for identifying significant breaks (default 0.05).
#' @param EFD Numeric scalar. Expected number of false discoveries. If provided, overrides \code{FDR}.
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
clust_BayesFuseMR <- function(fit, FDR = 0.05, EFD = NULL) {
  
  # Order theta by posterior means and define delta as their successive differences
  theta_samples <- fit$posterior$draws(variables = "theta", format = "draws_matrix") %>%
                      t()
  
  breaks <- t(apply(theta_samples[order(rowMeans(theta_samples)), ], 2, diff)) %>%
    as.data.frame() %>%
    setNames(paste0("delta[", seq_len(ncol(.)), "]")) %>%
    purrr::map_dfr(~{
      tibble(
        mean  = mean(.x),
        se    = sd(.x),
        pep   = 2 * min(mean(.x < 0), mean(.x > 0)),
        q25  = quantile(.x, 0.25),
        q75  = quantile(.x, 0.75)
      )
    }, .id = "variable") %>%
    arrange(pep) %>%
    mutate(qvalue = cummean(pep), qvalue2 = row_number() * qvalue) 
  
  if (!is.null(EFD)) {
    breaks <- filter(breaks, qvalue2 < EFD)
  } else {
    breaks <- filter(breaks, qvalue < FDR)
  }
  
  # Extract index of significant breaks from variable names
  break_indices <- breaks %>%
    mutate(index = as.integer(gsub("delta\\[(\\d+)\\]", "\\1", variable))) %>%
    arrange(index) %>%
    pull(index)
  
  # Assign cluster IDs
  theta_df <- fit$posterior$summary(variables = "theta", mean, sd) %>%
    mutate(snp = fit$snps) %>%
    mutate(index = as.integer(gsub("theta\\[(\\d+)\\]", "\\1", variable))) %>%
    slice(order(rowMeans(theta_samples))) %>%
    rename(theta_mean = mean, theta_se = sd) %>%
    mutate(cluster = cut(
      row_number(),
      breaks = c(0, break_indices, n()),  # breakpoints define cluster boundaries
      labels = FALSE,
      right = TRUE
    )) %>%
    select(snp, theta_mean, theta_se, cluster)
  
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
    dplyr::select(snp, cluster, cluster_size, cluster_mean, cluster_se, cluster_lcl, cluster_ucl, theta_mean, theta_se)
  
  return(theta_summary)
}
