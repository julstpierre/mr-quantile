#' Dummy dataset for MR examples
#'
#' A small simulated dataset used to illustrate the quantile Mendelian
#' randomization estimator implemented in the package.
#'
#' @format A data frame with 30 rows and 4 variables:
#' \describe{
#'   \item{beta_hat_y}{Estimated SNP-outcome association}
#'   \item{beta_hat_x}{Estimated SNP-exposure association}
#'   \item{sigma_y}{Standard error of the SNP-outcome estimate}
#'   \item{sigma_x}{Standard error of the SNP-exposure estimate}
#' }
#'
#' @source Simulated data
"dummy_data"