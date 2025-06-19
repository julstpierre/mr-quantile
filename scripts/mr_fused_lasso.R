library(CVXR)
library(genlasso)
library(dplyr)

mr_fused_lasso <- function(beta_hat_y, beta_hat_x , sigma2_y, nlambda){
  
  n <- length(beta_hat_y)
  
  # Calculate ratios and weights, and reorder IVs
  r <- beta_hat_y / beta_hat_x
  var_r <- sigma2_y / beta_hat_x^2
  w = 1 / sqrt(var_r)
  r_w <- (w * r)[order(r)]
  X_w <- cbind(diag(w[order(r)]), diag((w / beta_hat_x)[order(r)]))
  
  # Use genlasso package to obtain grid for lambda1
  genlasso.fit <- genlasso(r_w, X = X_w[, 1:n], D = diff(diag(n)))
  lambda1_seq <- genlasso.fit$lambda
  
  # Calculate grid for lambda2
  lambda2_max <- apply(abs(X_w[,(n+1):(2*n)] %*% (r_w - genlasso.fit$fit)), 2, max)
  lambda2_seq <- lapply(1:length(lambda2_max), function(i) exp(seq(log(lambda2_max[i]), log(0.01*lambda2_max[i]), length.out = nlambda))) %>% unlist()
  
  # Combine lambda1 and lambda 2 into same grid and fit fused lasso
  lambda_grid <- cbind(lambda1 = rep(lambda1_seq, each = nlambda), lambda2 = lambda2_seq)
  lasso.fit <- apply(lambda_grid, 1, function(grid) fused_lasso_fit(y = r_w, X = X_w, lambda1 = grid["lambda1"], lambda2 = grid["lambda2"]))

  list(theta_hat = sapply(lasso.fit, function(fit) fit$theta_hat),
       alpha_hat = sapply(lasso.fit, function(fit) fit$alpha_hat),
       fit = sapply(lasso.fit, function(fit) fit$fit),
       y = r_w,
       lambda_grid = lambda_grid
       )
  
}

fused_lasso_fit <- function(y, X, lambda1, lambda2, ZERO_THRESHOLD = 1e-6){
  
  n <- length(y)
  
  ## Define_Parameters
  theta <- Variable(n)
  alpha <- Variable(n)
  
  ## Define penalties
  fused_penalty <- p_norm(diff(theta), 1) #lambda1
  lasso_penalty <- p_norm(alpha, 1) #lambda2
  
  # Define regression model
  y_hat <- X[, 1:n] %*% theta + X[, (n+1):(2*n)] %*% alpha
  
  # Objective function
  objective <- Minimize(0.5*sum_squares(y - y_hat) + lambda1 * fused_penalty + lambda2 * lasso_penalty)
  fit <- solve(Problem(objective))
  
  theta_hat <- fit$getValue(theta)
  alpha_hat <- fit$getValue(alpha)
  
  ## Zero out stuff before returning
  theta_hat[abs(theta_hat) < ZERO_THRESHOLD] <- 0.0
  alpha_hat[abs(alpha_hat) < ZERO_THRESHOLD] <- 0.0
  list(theta_hat = theta_hat,
       alpha_hat = alpha_hat,
       fit = X[, 1:n] %*% theta_hat + X[, (n+1):(2*n)] %*% alpha_hat
       )
  
}
