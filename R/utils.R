#' Quantile regression check function
#'
#' Computes the check (pinball) loss used in quantile regression.
#'
#' @param u Numeric vector of residuals.
#' @param tau Quantile level.
#'
#' @return Numeric vector of loss values.
#'
#' @keywords internal
rho_tau <- function(u, tau) {
  u * (tau - (u < 0))
}


#' Weighted quantile
#'
#' Computes a weighted quantile using linear interpolation.
#'
#' @param x Numeric vector.
#' @param w Non-negative weights.
#' @param prob Quantile level.
#'
#' @return Weighted quantile of `x`.
#'
#' @keywords internal
weighted_quantile <- function(x, w, prob = 0.5) {
  
  o <- order(x)
  x <- x[o]
  w <- w[o]
  
  cw <- cumsum(w) - 0.5 * w
  cw <- cw / sum(w)
  
  idx <- which(cw < prob)
  
  if (length(idx) == 0) return(x[1])
  if (length(idx) == length(x)) return(x[length(x)])
  
  k <- max(idx)
  
  ratio <- (prob - cw[k]) / (cw[k + 1] - cw[k])
  x[k] + (x[k + 1] - x[k]) * ratio
}
