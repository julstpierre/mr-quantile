#--------------------------------------------------------------------------------
# Simulate data
#--------------------------------------------------------------------------------
sim_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(theta_1 = 0.2, theta_2 = 0, theta_3 = -1, nclust = 4, p = 50, h2 = 0.10, n = 50000)
  n <- args$n
  theta_1 <- args$theta_1
  theta_3 <- args$theta_3
  nclust <- args$nclust
  p <- args$p
  h2 <- args$h2
  
  if (nclust == 1){
    theta_2 <- rep(args$theta_2, p)
  } else if (nclust == 4){
      if (p == 30){
        theta_2 <- rep(c(args$theta_2, -0.4, 0, 0.4), times = c(12, 6, 6, 6))
      } else if (p == 50){
        theta_2 <- rep(c(args$theta_2, -0.4, 0, 0.4), times = c(20, 12, 10, 8))
      } else if (p == 100){
        theta_2 <- rep(c(args$theta_2, -0.4, 0, 0.4), times = c(45, 20, 20, 15))
      }
  }
  
  maf <- runif(p, 0.05, 0.5)
  beta_x_raw <- rnorm(p, mean = 0, sd = 1)
  beta_x <- sqrt( h2 / sum( beta_x_raw^2 * 2 * maf * ( 1 - maf ) ) ) * beta_x_raw
  beta_m <- theta_1 * beta_x + rnorm(p, sd = 1)
  beta_y <- theta_2 * beta_x + theta_3 * beta_m
  
  sigma2 <-  1 / ( 2 * n * maf * ( 1 - maf ) ) 
  beta_hat_x <- rnorm( p, beta_x, sqrt( sigma2 ) )
  beta_hat_m <- rnorm( p, beta_m, sqrt( sigma2 ) )
  beta_hat_y <- rnorm( p, beta_y, sqrt( sigma2 ) )
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2, 
              "beta_hat_m" = beta_hat_m, 
              "sigma2_m" = sigma2,
              "theta_T" = theta_2 + theta_1 * theta_3,
              "theta_D" = theta_2,
              "nclust" = nclust,
              "p" = p,
              "maf" = maf
              )
         )
}

sim_indv_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(theta_1 = 0.2, theta_2 = 0, theta_3 = -1, nclust = 4, p = 50, h2 = 0.10, n = 50000)
  n <- args$n
  theta_1 <- args$theta_1
  theta_3 <- args$theta_3
  nclust <- args$nclust
  p <- args$p
  h2 <- args$h2
  
  if (nclust == 1){
    theta_2 <- rep(args$theta_2, p)
  } else {
    if (p == 30){
      theta_2 <- rep(c(args$theta_2, -0.4, 0.8 - args$theta_2, 0.4), times = c(12, 6, 6, 6))
    } else if (p == 50){
      theta_2 <- rep(c(args$theta_2, -0.4, 0.8 - args$theta_2, 0.4), times = c(20, 10, 10, 10))
    } else if (p == 100){
      theta_2 <- rep(c(args$theta_2, -0.4, 0.8 - args$theta_2, 0.4), times = c(45, 20, 20, 15))
    }
  }
  
  # Simulate genetic data and confounders
  maf <- runif(p, 0.05, 0.5)
  G <- sapply(maf, function(maf) rbinom(n, 2, maf))
  U <- rnorm(n)
  epsilon_w <- mvrnorm(n, mu = rep(0, nclust), Sigma = diag(1, nclust))
  epsilon_x <- rnorm(n)
  epsilon_m <- rnorm(n)
  epsilon_y <- rnorm(n)
  
  # Simulate latent factors with 20% of variance explained by genetic predictors
  delta_w <- mvrnorm(p, mu = rep(0.2, nclust), Sigma = diag(0.1^2, nclust)) * sapply(1:nclust, function(i) theta_2 == unique(theta_2)[i])
  h2_w <- 0.2
  scaling_factors <- sqrt( h2_w / (1 - h2_w) * apply(epsilon_w, 2, var) / apply(G %*% delta_w, 2, var))
  if (nclust == 1){
    delta_w <- delta_w * scaling_factors
  } else {
    delta_w <- delta_w %*% diag(scaling_factors)
  }
  W <- G %*% delta_w + epsilon_w
  apply(G %*% delta_w, 2, var) / apply(W, 2, var) # ~0.2
  
  # Simulate exposure
  delta_x <- runif(nclust, -1, 1)
  phi <- h2 * ( var( U ) + var( epsilon_x ) ) / ( var( G %*% ( delta_w %*% delta_x ) ) - h2 * var( W %*% delta_x ) )
  delta_x <- delta_x %*% sqrt(phi)
  X <- W %*% delta_x + U + epsilon_x
  var(G %*% delta_w %*% delta_x) / var(X) # ~0.1
  
  # Simulate mediator
  delta_m <- rnorm(p, mean = 0, sd = 0.1)
  M <- theta_1 * X + U + epsilon_m + G %*% delta_m
  
  # Simulate outcome
  gamma <- delta_w %*% delta_x
  alpha <- gamma * (theta_2 - theta_2[1])
  delta_y <- qr.solve(delta_w, alpha)
  Y <- theta_2[1] * X + W %*% delta_y + theta_3 * M + U + epsilon_y 
  
  # Summary stats for X
  sampX = seq(1, n)
  Gc <- scale(G[sampX, ], scale = FALSE)
  Xc <- X[sampX] - mean(X[sampX])
  XtG <- crossprod(Gc, Xc)            
  GtG <- colSums(Gc^2)                   
  beta_hat_x <- XtG / GtG
  residuals2 <- colSums((matrix(Xc, n, p) - Gc * rep(beta_hat_x, each = n))^2)
  sigma2_x <- residuals2 / (n - 2) / GtG
  
  # Summary stats for X
  sampM = seq(1, n)
  Gc <- scale(G[sampM, ], scale = FALSE)
  Mc <- M[sampM] - mean(M[sampM])
  XtG <- crossprod(Gc, Mc)            
  GtG <- colSums(Gc^2)                   
  beta_hat_m <- XtG / GtG
  residuals2 <- colSums((matrix(Mc, n, p) - Gc * rep(beta_hat_m, each = n))^2)
  sigma2_m <- residuals2 / (n - 2) / GtG
  
  # Summary stats for Y
  sampY <- seq(1, n)
  Yc <- Y[sampY] - mean(Y[sampY])
  Gc <- scale(G[sampY, ], scale = FALSE)
  XtG <- crossprod(Gc, Yc)
  GtG <- colSums(Gc^2)
  beta_hat_y <- XtG / GtG
  residuals2 <- colSums((matrix(Yc, n, p) - Gc * rep(beta_hat_y, each = n))^2)
  sigma2_y <- residuals2 / (n - 2) / GtG
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2_x, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2_y, 
              "beta_hat_m" = beta_hat_m, 
              "sigma2_m" = sigma2_m,
              "theta_T" = theta_2 + theta_1 * theta_3,
              "theta_D" = theta_2,
              "nclust" = nclust,
              "p" = p,
              "maf" = maf
  )
  )
}

#--------------------------------------------------------------------------------
# Fit MR_IVW on simulated data
#--------------------------------------------------------------------------------
theta_ivw <- function(D, snps = 1:D$p){

  total.effect <- lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1)$coef[1]
  resid.std.error <-summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$sigma
  se.total.effect <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$coef[1,2]
  ci.upper.total <- max(total.effect + qnorm(0.975) * se.total.effect / resid.std.error,
                       total.effect + qt(0.975, df=length(snps)-1) * se.total.effect)
  ci.lower.total <- min(total.effect - qnorm(0.975) * se.total.effect / resid.std.error,
                       total.effect - qt(0.975, df=length(snps)-1) * se.total.effect)

  return(cbind("Estimate" = total.effect, "StdError" = se.total.effect, "LCL" = ci.lower.total, "UCL" = ci.upper.total))
}

theta_mvivw <- function(D, snps = 1:D$p){

  direct.effect <- lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] + D$beta_hat_m[snps] - 1, weights = D$sigma2_y[snps]^-1)$coef[1]
  resid.std.error <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$sigma
  se.direct.effect <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] + D$beta_hat_m[snps] - 1, weights = D$sigma2_y[snps]^-1))$coef[1,2]
  ci.upper.direct <- max(direct.effect + qnorm(0.975) * se.direct.effect / resid.std.error,
                       direct.effect + qt(0.975, df=length(snps)-1) * se.direct.effect)
  ci.lower.direct <- min(direct.effect - qnorm(0.975) * se.direct.effect / resid.std.error,
                       direct.effect - qt(0.975, df=length(snps)-1) * se.direct.effect)

  return(cbind("Estimate" = direct.effect, "StdError" = se.direct.effect, "LCL" = ci.lower.direct, "UCL" = ci.upper.direct))
}

fit_ivw <- function(D, snps = 1:D$p){
  theta_T <- theta_ivw(D, snps)
  theta_D <- theta_mvivw(D, snps)
  theta <- cbind(
          theta_T_hat = theta_T[,"Estimate"], theta_T_se = theta_T[,"StdError"],
          theta_T_lcl = theta_T[,"LCL"], theta_T_ucl = theta_T[,"UCL"],
          theta_D_hat = theta_D[,"Estimate"], theta_D_se = theta_D[,"StdError"],
          theta_D_lcl = theta_D[,"LCL"], theta_D_ucl = theta_D[,"UCL"]
        )
  row.names(theta) <- NULL
  return(theta)
}

# theta_ivw <- function(D, snps = 1:D$p){
#   MRInputObject = mr_input(bx = D$beta_hat_x[snps], bxse = sqrt(D$sigma2_x[snps]), by = D$beta_hat_y[snps], byse = sqrt(D$sigma2_y[snps]))
#   MR_IVW <- mr_ivw(MRInputObject)
#   return(cbind("Estimate" = unname(MR_IVW$Estimate[1]), "StdError" = unname(MR_IVW$StdError[1]), "Pvalue" = unname(MR_IVW$Pvalue[1])))
# }
# 
# theta_mvivw <- function(D, snps = 1:D$p){
#   MRInputObject = mr_mvinput(bx = cbind(D$beta_hat_x[snps], D$beta_hat_m[snps]), 
#                              bxse = sqrt(cbind(D$sigma2_x[snps], D$sigma2_m[snps])), 
#                              by = D$beta_hat_y[snps], 
#                              byse = sqrt(D$sigma2_y[snps])
#   )
#   MR_IVW <- mr_mvivw(MRInputObject)
#   return(cbind("Estimate" = unname(MR_IVW$Estimate[1]), "StdError" = unname(MR_IVW$StdError[1]), "Pvalue" = unname(MR_IVW$Pvalue[1])))
# }
# 
# fit_ivw <- function(D, snps = 1:D$p){
#   theta_T <- theta_ivw(D, snps)
#   theta_D <- theta_mvivw(D, snps)
#   theta <- cbind(
#     theta_T_hat = theta_T[,"Estimate"], theta_T_se = theta_T[,"StdError"], 
#     theta_T_lcl = theta_T[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
#     theta_T_ucl = theta_T[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
#     theta_D_hat = theta_D[,"Estimate"], theta_D_se = theta_D[,"StdError"],
#     theta_D_lcl = theta_D[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"],
#     theta_D_ucl = theta_D[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"]
#   )
#   row.names(theta) <- NULL
#   return(theta)
# }

#--------------------------------------------------------------------------------
# Fit MR_Egger on simulated data
#--------------------------------------------------------------------------------
theta_egger <- function(D, snps = 1:D$p){
    MRInputObject = mr_input(bx = D$beta_hat_x[snps], bxse = sqrt(D$sigma2_x[snps]), by = D$beta_hat_y[snps], byse = sqrt(D$sigma2_y[snps]))
    MR_Egger <- mr_egger(MRInputObject)
    return(cbind("Estimate" = unname(MR_Egger$Estimate[1]), "StdError" = unname(MR_Egger$StdError.Est[1]), "Pvalue" = unname(MR_Egger$Pvalue.Est[1])))
}

theta_mvegger <- function(D, snps = 1:D$p){
  MRInputObject = mr_mvinput(bx = cbind(D$beta_hat_x[snps], D$beta_hat_m[snps]), 
                             bxse = sqrt(cbind(D$sigma2_x[snps], D$sigma2_m[snps])), 
                             by = D$beta_hat_y[snps], 
                             byse = sqrt(D$sigma2_y[snps])
  )
  MR_Egger <- mr_mvegger(MRInputObject)
  return(cbind("Estimate" = unname(MR_Egger$Estimate[1]), "StdError" = unname(MR_Egger$StdError.Est[1]), "Pvalue" = unname(MR_Egger$Pvalue.Est[1])))
}

fit_egger <- function(D, snps = 1:D$p){
  theta_T <- theta_egger(D, snps)
  theta_D <- theta_mvegger(D, snps)
  theta <- cbind(
    theta_T_hat = theta_T[,"Estimate"], theta_T_se = theta_T[,"StdError"], 
    theta_T_lcl = theta_T[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
    theta_T_ucl = theta_T[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
    theta_D_hat = theta_D[,"Estimate"], theta_D_se = theta_D[,"StdError"],
    theta_D_lcl = theta_D[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"],
    theta_D_ucl = theta_D[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"]
  )
  row.names(theta) <- NULL
  return(theta)
}

#--------------------------------------------------------------------------------
# Fit MR_Median on simulated data
#--------------------------------------------------------------------------------
theta_median <- function(D, snps = 1:D$p){
    MR_Median <- summary(rq(D$beta_hat_y[snps]~-1+D$beta_hat_x[snps], weights = D$sigma2_y[snps]^-1), se = "boot")$coef
    return(cbind("Estimate" = MR_Median[1], "StdError" = MR_Median[2], "Pvalue" = MR_Median[4]))
}

theta_mvmedian <- function(D, snps = 1:D$p){
  MR_Median <- summary(rq(D$beta_hat_y[snps]~-1+D$beta_hat_x[snps] + D$beta_hat_m[snps], weights = D$sigma2_y[snps]^-1), se = "boot")$coef[1, ]
  return(cbind("Estimate" = unname(MR_Median[1]), "StdError" = unname(MR_Median[2]), "Pvalue" = unname(MR_Median[4])))
}

fit_median <- function(D, snps = 1:D$p){
  theta_T <- theta_median(D, snps)
  theta_D <- theta_mvmedian(D, snps)
  theta <- cbind(
    theta_T_hat = theta_T[,"Estimate"], theta_T_se = theta_T[,"StdError"], 
    theta_T_lcl = theta_T[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
    theta_T_ucl = theta_T[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_T[,"StdError"],
    theta_D_hat = theta_D[,"Estimate"], theta_D_se = theta_D[,"StdError"],
    theta_D_lcl = theta_D[,"Estimate"] - qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"],
    theta_D_ucl = theta_D[,"Estimate"] + qt(0.975, df = length(snps) - 1) * theta_D[,"StdError"]
  )
  return(theta)
}

#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
fit_mrclust <- function(D){
    ratios <- D$beta_hat_y / D$beta_hat_x
    ratios.se <- sqrt(D$sigma2_y) / abs(D$beta_hat_x)
    res_em <- mr_clust_em(theta = ratios, 
                          theta_se = ratios.se, 
                          bx = D$beta_hat_x, 
                          by = D$beta_hat_y, 
                          bxse = sqrt(D$sigma2_x), 
                          byse = sqrt(D$sigma2_y),
                          obs_names = paste0("SNP", 1:D$p)
                         )
    
    res_em_best <- res_em$results$best
    table(res_em_best$cluster_mean)
    
    # Find cluster with the most number of SNPs
    main_cluster <- names(which.max(table(res_em_best$cluster_mean)))
    snps <- filter(res_em_best, cluster_mean == main_cluster) %>%
            mutate(rows = substr(observation, 4, 6)) %>%
            pull(rows) %>%
            as.numeric()
    
    # Fit IVW or median estimator on main cluster
    nclust <- length(unique(res_em_best$cluster))
    if (nclust == 1) {
      theta_main <- theta_ivw(D)
    } else {
      theta_main <- theta_median(D, snps)
    }
    
    return(
      cbind(
           theta_main,
           nclust = nclust,
           adj.rand = adjustedRandIndex(res_em_best$cluster, D$theta_T),
           rand = rand.index(res_em_best$cluster, D$theta_T)
           )
      )
}

fit_mvmrclust <- function(D){
  fit <- lm(D$beta_hat_y ~ D$beta_hat_m, weights = D$sigma2_y^-1)
  resid.y <- residuals(fit)
  resid.se <- sqrt(D$sigma2_y * (1 - hatvalues(fit)))
  ratios <- resid.y / D$beta_hat_x
  ratios.se <- resid.se / abs(D$beta_hat_x)
  res_em <- mr_clust_em(theta = ratios, 
                        theta_se = ratios.se, 
                        bx = D$beta_hat_x, 
                        by = resid.y, 
                        bxse = sqrt(D$sigma2_x), 
                        byse = resid.se,
                        obs_names = paste0("SNP", 1:D$p)
  )
  
  res_em_best <- res_em$results$best
  table(res_em_best$cluster_mean)
  
  # Find cluster with the most number of SNPs
  main_cluster <- names(which.max(table(res_em_best$cluster_mean)))
  snps <- filter(res_em_best, cluster_mean == main_cluster) %>%
    mutate(rows = substr(observation, 4, 6)) %>%
    pull(rows) %>%
    as.numeric()
  
  # Fit IVW or median estimator on main cluster
  nclust <- length(unique(res_em_best$cluster))
  if (nclust == 1) {
    theta_main <- fit_ivw(D)
  } else {
    theta_main <- fit_median(D, snps)
  }
  
  return(
    cbind(
      theta_main,
      nclust = nclust,
      adj.rand = adjustedRandIndex(res_em_best$cluster, D$theta_T),
      rand = rand.index(res_em_best$cluster, D$theta_T)
    )
  )
}

#--------------------------------------------------------------------------------
# Fit generalized fused lasso on simulated data
#--------------------------------------------------------------------------------
fit_fused_lasso <- function(D, k = 1, gamma = 2, ebic_coef = 7, post_hoc = "mean", lambda_list = NULL, max_clusters = NULL){
  
  n <- length(D$beta_hat_y)
  if (is.null(max_clusters)) { max_clusters <- n/2 }
  y <- D$beta_hat_y/sqrt(D$sigma2_y)
  X <- diag(D$beta_hat_x/sqrt(D$sigma2_y))
  
  # For each pair (i, j), we want a row with 1 at position i, -1 at position j
  pairs <- do.call(rbind, lapply(1:(n - 1), function(i) cbind(i, (i + 1):min(n, i + k))))
  m <- nrow(pairs)
  Dmat <- sparseMatrix(
    i = rep(seq_len(m), 2),
    j = c(pairs[,1], pairs[,2]),
    x = c(rep(1, m), rep(-1, m)),
    dims = c(m, n)
  ) %>% as.matrix()
  
  # Apply a l2 regularization to the ratios
  if (is.null(lambda_list)) { lambda_list <- seq(0, 3, 0.1) }
  results <- vector("list", length(lambda_list))
  
  for (i in seq_along(lambda_list)) {
      lambda <- lambda_list[i]
      ytilde <- c(y, rep(0, n))
      Xtilde <- rbind(X, sqrt(lambda) * diag(n))
      ridge.ratios <- lm(ytilde ~ -1 + Xtilde)$coefficients %>% as.vector()
      
      # Order data by regularized ratios
      ord <- order(ridge.ratios)
      ytilde <- c(y[ord], rep(0, n))
      Xtilde <- rbind(X[ord, ord], sqrt(lambda) * diag(n))
      
      # Compute adaptive weights using differences in ratios
      eps     <- 0    # stability constant
      diffs   <- sweep(Dmat, 2, ridge.ratios[ord], `*`) %>% rowSums() %>% abs()
      weights <- 1 / (diffs^gamma + eps)
      
      # Fit model
      genlasso.fit <- genlasso(ytilde, Xtilde, D = diag(weights) %*% Dmat)
      
      # Compute EBIC for each lambda
      rss <- colSums((genlasso.fit$fit[1:n, ] - ytilde[1:n])^2)
      df <- genlasso.fit$df  # effective degrees of freedom
      ebic <- rss / n + log(n) * df / n + 2 * ebic_coef * log(choose(n, df)) / n
        
      # Return best fit
      results[[i]] <- list(
        lambda = lambda,
        theta_hat = round(genlasso.fit$beta[, which.min(ebic[df < max_clusters])], 4),
        ebic = min(ebic[df < max_clusters]),
        ord = ord,
        ratios = ridge.ratios
      )
  }
  
  # Select best model
  best_fit <- results[[which.min(sapply(results, function(x) x$ebic))]]
  best_fit$lambda
  table(best_fit$theta_hat)
  
  # Assign clusters based on distance criteria
  if (!is.null(post_hoc)){
    
    # Compute distance matrix
    cluster_means <- sort(unique(best_fit$theta_hat))
    names(cluster_means) <- cluster_means
    ratios.se <- sqrt(D$sigma2_y) / abs(D$beta_hat_x)
    distance_list <- lapply(results, function(res) {
      sapply(cluster_means, function(mu) {
        (res$ratios - mu)^2 / ratios.se^2
      })
    })
    
    # Assign clusters based on minimum distance
    if (post_hoc == "mean") {
      cluster <- apply(apply(simplify2array(distance_list), c(1, 2), mean), 1, which.min)
    } else if (post_hoc == "mode") {
      cluster <- apply(apply(simplify2array(distance_list), c(1, 3), which.min), 1, get.mode)
    }
    names(cluster) <- cluster_means[cluster]
    theta_hat <- as.numeric(cluster_means[cluster])
    table(theta_hat)
    
  } else {
    # Assign clusters based on fusion
    cluster_means <- sort(unique(best_fit$theta_hat))
    cluster <- match(best_fit$theta_hat[order(best_fit$ord)], cluster_means)
    names(cluster) <- cluster_means[cluster]
    theta_hat <- as.numeric(cluster_means[cluster])
    table(theta_hat)
  }
  
  # Find cluster with the most number of SNPs
  main_cluster <- table(cluster) |> which.max() |> names() |> as.numeric()
  snps <- which(cluster == main_cluster)
  
  # Fit IVW or median estimator on main cluster
  nclust <- length(unique(cluster))
  if (nclust == 1) {
    theta_main <- theta_ivw(D)
  } else {
    theta_main <- theta_median(D, snps)
  }
  
  return(
    cbind(
      theta_main,
      nclust = nclust,
      adj.rand = adjustedRandIndex(theta_hat, D$theta_T),
      rand = rand.index(theta_hat, D$theta_T)
    )
  )
  
}

fit_mvfused_lasso <- function(D, k = 1, gamma = 2, ebic_coef = 8, post_hoc = "mean", lambda_list = NULL, max_clusters = NULL){
  
  n <- length(D$beta_hat_y)
  if (is.null(max_clusters)) { max_clusters <- n/2 }
  y <- D$beta_hat_y/sqrt(D$sigma2_y)
  X <- diag(D$beta_hat_x/sqrt(D$sigma2_y))
  M <- D$beta_hat_m/sqrt(D$sigma2_y)
    
  # For each pair (i, j), we want a row with 1 at position i, -1 at position j
  pairs <- do.call(rbind, lapply(1:(n - 1), function(i) cbind(i, (i + 1):min(n, i + k))))
  m <- nrow(pairs)
  Dmat <- sparseMatrix(
    i = rep(seq_len(m), 2),
    j = c(pairs[,1], pairs[,2]),
    x = c(rep(1, m), rep(-1, m)),
    dims = c(m, n)
  ) %>% as.matrix()
  
  # Apply a l2 regularization to the ratios
  if (is.null(lambda_list)) { lambda_list <- seq(0, 4, 0.1) }
  results <- vector("list", length(lambda_list))
  
  for (i in seq_along(lambda_list)) {
    lambda <- lambda_list[i]
    ytilde <- c(y, rep(0, n + 1))
    Xtilde <- rbind(cbind(X, M), sqrt(lambda) * diag(n + 1))
    ridge.ratios <- lm(ytilde ~ -1 + Xtilde)$coefficients[1:n] %>% as.vector()
    
    # Order data by regularized ratios
    ord <- order(ridge.ratios)
    ytilde <- c(y[ord], rep(0, n + 1))
    Xtilde <- rbind(cbind(M[ord], X[ord, ord]), sqrt(lambda) * diag(n + 1))
      
    # Compute adaptive weights using differences in ratios
    eps     <- 0    # stability constant
    diffs   <- sweep(Dmat, 2, ridge.ratios[ord], `*`) %>% rowSums() %>% abs()
    weights <- 1 / (diffs^gamma + eps)
    
    # Fit model
    genlasso.fit <- genlasso(ytilde, Xtilde, D = cbind(0, diag(weights) %*% Dmat))
    
    # Compute EBIC for each lambda
    rss <- colSums((genlasso.fit$fit[1:n, ] - ytilde[1:n])^2)
    df <- genlasso.fit$df  # effective degrees of freedom
    ebic <- rss / n + log(n) * df / n + 2 * ebic_coef * log(choose(n, df)) / n
    
    # Return best fit
    results[[i]] <- list(
      lambda = lambda,
      theta_hat = round(genlasso.fit$beta[-1, which.min(ebic[df < max_clusters])], 4),
      ebic = min(ebic[df < max_clusters]),
      ord = ord,
      ratios = ridge.ratios, 
      sigma = unname(sqrt((rss / (n - df))[which.min(ebic[df < max_clusters])]))
    )
  }
  
  # Select best model
  best_fit <- results[[which.min(sapply(results, function(x) x$ebic))]]
  best_fit$lambda
  table(best_fit$theta_hat)
  
  # Assign clusters based on distance criteria
  if (!is.null(post_hoc)){
    
    # Compute distance matrix
    cluster_means <- sort(unique(best_fit$theta_hat))
    names(cluster_means) <- cluster_means
    ratios.se <- sqrt(D$sigma2_y) / abs(D$beta_hat_x)
    distance_list <- lapply(results, function(res) {
      sapply(cluster_means, function(mu) {
        (res$ratios - mu)^2 / ratios.se^2
      })
    })
    
    # Assign clusters based on minimum distance
    if (post_hoc == "mean") {
      cluster <- apply(apply(simplify2array(distance_list), c(1, 2), mean), 1, which.min)
    } else if (post_hoc == "mode") {
      cluster <- apply(apply(simplify2array(distance_list), c(1, 3), which.min), 1, get.mode)
    }
    names(cluster) <- cluster_means[cluster]
    theta_hat <- as.numeric(cluster_means[cluster])
    table(theta_hat)
    
  } else {
    # Assign clusters based on fusion
    cluster_means <- sort(unique(best_fit$theta_hat))
    cluster <- match(best_fit$theta_hat[order(best_fit$ord)], cluster_means)
    names(cluster) <- cluster_means[cluster]
    theta_hat <- as.numeric(cluster_means[cluster])
    table(theta_hat)
  }
  
  # Find cluster with the most number of SNPs
  main_cluster <- table(cluster) |> which.max() |> names() |> as.numeric()
  snps <- which(cluster == main_cluster)
  
  # Fit IVW or median estimator on main cluster
  nclust <- length(unique(cluster))
  if (nclust == 1) {
    theta_main <- fit_ivw(D)
  } else {
    theta_main <- fit_median(D, snps)
  }
  
  return(
    cbind(
      theta_main,
      nclust = nclust,
      adj.rand = adjustedRandIndex(theta_hat, D$theta_T),
      rand = rand.index(theta_hat, D$theta_T)
    )
  )
  
}
