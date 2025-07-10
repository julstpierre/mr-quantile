#--------------------------------------------------------------------------------
# Simulate data
#--------------------------------------------------------------------------------
sim_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(beta = 0.4, nclust = 4, p = 30, h2 = 0.30)
  n <- 50000
  beta <- args$beta
  nclust <- args$nclust
  p <- args$p
  h2 <- args$h2
  
  if (nclust == 1){
    theta <- rep(beta, p)
  } else if (nclust == 4){
      if (p == 30){
        theta <- rep(c(beta, -0.4, 0.8, 0), times = c(12, 6, 6, 6))
      } else if (p == 100){
        theta <- rep(c(beta, -0.4, 0.8, 0), times = c(45, 20, 20, 15))
      }
  }
  
  maf <- runif(p, 0.05, 0.5)
  beta_x <- rnorm( p, 0, sqrt( h2 / ( 2 * p * maf * (1 - maf) ) ) )
  beta_y <- theta * beta_x
  sigma2 <-  1 / ( 2 * n * maf * ( 1 - maf ) ) 
  
  beta_hat_x <- rnorm( p, beta_x, sqrt( sigma2 ) )
  beta_hat_y <- rnorm( p, beta_y, sqrt( sigma2 ) )
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2, 
              "beta" = beta,
              "theta" = theta,
              "nclust" = nclust,
              "p" = p
              )
         )
}

#--------------------------------------------------------------------------------
# Fit MR_IVW on simulated data
#--------------------------------------------------------------------------------
theta_ivw <- function(D, snps = 1:D$p){
    MRInputObject = mr_input(bx = D$beta_hat_x[snps], bxse = sqrt(D$sigma2_x[snps]), by = D$beta_hat_y[snps], byse = sqrt(D$sigma2_y[snps]))
    MR_IVW <- mr_ivw(MRInputObject)
    return(cbind("Estimate" = unname(MR_IVW$Estimate[1]), "StdError" = unname(MR_IVW$StdError[1]), "Pvalue" = unname(MR_IVW$Pvalue[1])))
}

#--------------------------------------------------------------------------------
# Fit MR_Egger on simulated data
#--------------------------------------------------------------------------------
theta_egger <- function(D){
    MRInputObject = mr_input(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), by = D$beta_hat_y, byse = sqrt(D$sigma2_y))
    MR_Egger <- mr_egger(MRInputObject)
    return(cbind("Estimate" = unname(MR_Egger$Estimate[1]), "StdError" = unname(MR_Egger$StdError.Est[1]), "Pvalue" = unname(MR_Egger$Pvalue.Est[1])))
}

#--------------------------------------------------------------------------------
# Fit MR_Median on simulated data
#--------------------------------------------------------------------------------
theta_median <- function(D){
    p <- length(D$beta_hat_y)
    theta_Median <- rq(D$beta_hat_y~-1+D$beta_hat_x, weights = D$sigma2_y^-1)$coefficients[1]
    names(theta_Median)[1] <- "Estimate"
    theta_Median["Std. Error"] <- replicate(1000, {
      beta_hat_y_k <- rnorm(p, D$beta_hat_y, sqrt(D$sigma2_y))
      beta_hat_x_k <- rnorm(p, D$beta_hat_x, sqrt(D$sigma2_x))
      rq(beta_hat_y_k~-1+beta_hat_x_k, weights = D$sigma2_y^-1)$coefficients[1]
    }) %>% sd()
    theta_Median["Pr(>|t|)"] <- min(1, 2*pt(q=theta_Median["Estimate"] / theta_Median["Std. Error"], df = p-1, lower.tail = FALSE))
    return(theta_Median)
}

#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
theta_mrclust <- function(D){
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
    
    # Find cluster with the most number of SNPs
    main_cluster <- names(which.max(table(res_em_best$cluster)))
    snps <- filter(res_em_best, cluster == main_cluster) %>%
            mutate(rows = substr(observation, 4, 6)) %>%
            pull(rows) %>%
            as.numeric()

    # Fit IVW on main cluster
    return(
      cbind(
           theta_ivw(D, snps), 
           nclust = length(unique(res_em_best$cluster)),
           adj.rand = adjustedRandIndex(res_em_best$cluster, D$theta),
           rand = rand.index(res_em_best$cluster, D$theta)
           )
      )
}

#--------------------------------------------------------------------------------
# Fit generalized fused lasso on simulated data
#--------------------------------------------------------------------------------
theta_fused_lasso <- function(D, k=1){
  
  n <- length(D$beta_hat_y)
  ratios <- D$beta_hat_y / D$beta_hat_x
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
  lambda_list <- seq(0, 3, 0.1)
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
      gamma   <- 2     # can tune
      eps     <- 0    # stability constant
      diffs   <- abs(diff(ridge.ratios[ord]))
      weights <- 1 / (diffs^gamma + eps)
      
      # Fit model
      genlasso.fit <- genlasso(ytilde, Xtilde, D = diag(weights) %*% Dmat)
      
      # Compute EBIC for each lambda
      rss <- colSums((genlasso.fit$fit[1:n, ] - ytilde[1:n])^2)
      df <- genlasso.fit$df  # effective degrees of freedom
      ebic <- function(gamma) {
        rss / n + log(n) * df / n + 2 * 1 * log(choose(D$p-1, df-1)) / n
      }
      
      # Return best fit
      results[[i]] <- list(
        lambda = lambda,
        theta_hat = round(genlasso.fit$beta[, which.min(ebic(1))], 4),
        ebic = min(ebic(1)),
        ord = ord
      )
  }
  
  # Select best model
  best_fit <- results[[which.min(sapply(results, function(x) x$ebic))]]
  best_fit$lambda
  table(best_fit$theta_hat)
  
  # Find cluster with the most number of SNPs
  main_cluster <- names(which.max(table(best_fit$theta_hat)))
  snps <- ord[best_fit$theta_hat == main_cluster]
  
  # Fit IVW on main cluster
  return(
    cbind(
      theta_ivw(D, snps), 
      nclust = length(unique(best_fit$theta_hat)),
      adj.rand = adjustedRandIndex(best_fit$theta_hat, D$theta[best_fit$ord]),
      rand = rand.index(best_fit$theta_hat, D$theta[best_fit$ord])
    )
  )
  
}

#--------------------------------------------------------------------------------
# Fit Bayesian Fused lasso on simulated data
#--------------------------------------------------------------------------------
theta_BayesFuseMR <- function(D, FDR = NULL, EFD = NULL){
  
  fit <- fit_BayesFuseMR(
    beta_hat_y = D$beta_hat_y,
    beta_hat_x = D$beta_hat_x,
    sigma_y = sqrt(D$sigma2_y),
    sigma_x = sqrt(D$sigma2_x),
    snps = paste0("SNP", 1:D$p),
    seed = 123
  )
  
  results_list <- list()
    
  if (!is.null(FDR)) {
    # Loop over FDR values
    for (fdr in FDR) {
          best_fit <- clust_BayesFuseMR(fit, FDR = fdr)
            
          # Find cluster with the most number of SNPs
          main_cluster <- names(which.max(table(best_fit$cluster)))
          theta_main <- best_fit %>%
            distinct(cluster, .keep_all = TRUE) %>%
            filter(cluster == main_cluster) %>%
            dplyr::select(cluster_mean, cluster_se) %>%
            rename(Estimate = cluster_mean, StdError = cluster_se) %>%
            mutate(Pvalue = 2 * pnorm(abs(Estimate / StdError), lower.tail = FALSE))

          ord <- order(as.numeric(substr(fit$snps, 4, 6)))

          res <- cbind(
            theta_main, 
            nclust = length(unique(best_fit$cluster)),
            adj.rand = adjustedRandIndex(best_fit$cluster[ord], D$theta),
            rand = rand.index(best_fit$cluster[ord], D$theta),
            FDR = fdr
          )

          results_list[[as.character(fdr)]] <- res
        }
  } else if (!is.null(EFD)) {
    # Loop over EFD values
    for (efd in EFD) {
        best_fit <- clust_BayesFuseMR(fit, EFD = efd)

        # Find cluster with the most number of SNPs
        main_cluster <- names(which.max(table(best_fit$cluster)))
        theta_main <- best_fit %>%
            filter(cluster == main_cluster) %>%
            dplyr::select(cluster_mean, cluster_se) %>%
            rename(Estimate = cluster_mean, StdError = cluster_se) %>%
            mutate(Pvalue = 2 * pnorm(abs(Estimate / StdError), lower.tail = FALSE))

        ord <- order(as.numeric(substr(fit$snps, 4, 6)))

        res <- cbind(
            theta_main, 
            nclust = length(unique(best_fit$cluster)),
            adj.rand = adjustedRandIndex(best_fit$cluster[ord], D$theta),
            rand = rand.index(best_fit$cluster[ord], D$theta),
            EFD = efd
          )

          results_list[[as.character(efd)]] <- res
        }
  }
  
  # Combine all results into one data frame
  return(dplyr::bind_rows(results_list))

}
