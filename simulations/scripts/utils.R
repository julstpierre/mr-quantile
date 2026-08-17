#--------------------------------------------------------------------------------
# Simulate summary-level data
#--------------------------------------------------------------------------------

sim_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(p = 50, n = 20000, theta = 0.1, h2_x = 0.5, h2_y = 0.1, h2_u = 0.1, q = 0.6)
  p <- args$p
  n <- args$n
  theta <- args$theta
  h2_x <- args$h2_x
  h2_y <- args$h2_y
  h2_u <- args$h2_u
  q <- args$q
  m <- q * p
  invalid.IVs <- seq_len(m)
  
  # Simulate direct effects
  sigma2_x <- h2_x / p
  gamma <- rnorm(p, mean = 0, sd = sqrt(sigma2_x))
  
  # Simulate uncorrelated pleiotropic effects
  sigma2_alpha <- h2_y / p
  alpha <- rep(0, p)
  alpha[invalid.IVs] <- rnorm(m, mean = 0, sd = sqrt(sigma2_alpha))

  # Simulate correlated pleiotropic effects
  sigma2_phi <- h2_u / p
  phi <- rep(0, p)
  phi[invalid.IVs] <- rnorm(m, mean = 0, sd = sqrt(sigma2_phi))
  
  # Generate summary statistics
  sigma2 <-  rep(1 / n, p)
  beta_hat_x <- rnorm(p, gamma + phi, sqrt( sigma2 ))
  beta_hat_y <- rnorm(p, theta * ( gamma + phi ) + alpha + phi, sqrt( sigma2 ))
    
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2, 
              "theta" = theta,
              "invalid.IVs" = invalid.IVs,
              "n" = n, 
              "q" = q,
              "p" = p
              )
         )
}

#--------------------------------------------------------------------------------
# Simulate individual-level data
#--------------------------------------------------------------------------------
sim_indv_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(theta = 0.1, p = 30, n = 50000, q = 0.6, InSIDE = FALSE, pi = 0.7)
  n <- args$n 
  theta <- args$theta
  p <- args$p
  q <- args$q
  m <- q * p
  InSIDE <- args$InSIDE
  pi <- args$pi
  
  # Simulate genetic data, errors and invalid IVs
  maf <- runif(p, 0.1, 0.3)
  G <- sapply(maf, function(maf) rbinom(2 * n, 2, maf))
  epsilon_u <- rnorm(2 * n)
  epsilon_x <- rnorm(2 * n)
  epsilon_y <- rnorm(2 * n)
  invalid.IVs <- sample(1:p, m)
  
  # Simulate confounders
  phi <- rep(0, p)
  if (!InSIDE){ 
    phi[invalid.IVs] <- runif(m, -0.1, 0.1)
  }
  U <- G %*% phi + epsilon_u
  
  # Simulate exposure
  gamma <- runif(p, 0.1, 0.2) * sample(c(-1, 1), size = p, replace = TRUE, prob = c(pi, 1 - pi))
  X <- G %*% gamma + U + epsilon_x
  
  # F-statistic
  (R2 <- var(G %*% gamma) / var(X))
  R2 / (1 - R2) * (n - p - 1) / p
  
  # Simulate outcome
  alpha <- rep(0, p)
  alpha[invalid.IVs] <- runif(m, 0.2, 0.3)
  Y <- X %*% theta + G %*% alpha + U + epsilon_y
  
  # Summary stats for X
  sampX <- seq(1, n)
  regX <- apply(X[sampX, , drop = FALSE], 2, function(xcol) sumstat(xcol, G[sampX, ]))
  beta_hat_x <- as.vector(sapply(regX, `[[`, "beta_hat"))
  sigma2_x    <- as.vector(sapply(regX, `[[`, "sigma2"))
  
  # Summary stats for Y
  sampY <- seq(n+1, 2*n)
  regY <- sumstat(Y[sampY], G[sampY, ])
  beta_hat_y <- as.vector(regY$beta_hat)
  sigma2_y <- as.vector(regY$sigma2)
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2_x, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2_y,
              "theta" = theta,
              "n" = n,
              "p" = p,
              "q" = q,
              "invalid.IVs" = invalid.IVs,
              "alpha" = alpha,
              "phi" = phi,
              "gamma" = gamma
  )
  )
}

sumstat <- function(X, G){
  n <- nrow(G); p <- ncol(G)
  Xc <- X - mean(X)
  Gc <- scale(G, scale = FALSE)
  XtG <- crossprod(Gc, Xc)
  GtG <- colSums(Gc^2)
  beta_hat <- XtG / GtG
  residuals2 <- colSums((matrix(Xc, n, p) - Gc * rep(beta_hat, each = n))^2)
  sigma2 <- residuals2 / (n - 2) / GtG
  return(list(beta_hat = beta_hat, sigma2 = sigma2))
}

#--------------------------------------------------------------------------------
# Fit MR_IVW on simulated data
#--------------------------------------------------------------------------------
fit_mr_ivw <- function(D, snps = 1:D$p){
            TwoSampleMR::mr_ivw(D$beta_hat_x[snps], D$beta_hat_y[snps], 
                                sqrt(D$sigma2_x[snps]), sqrt(D$sigma2_y[snps])) %>% 
            as_tibble() %>% 
            mutate(Estimate = b,
                   StdError = se,
                   LCL = Estimate - qnorm(0.975) * StdError,
                   UCL = Estimate + qnorm(0.975) * StdError
            ) %>% 
            select(Estimate, StdError, LCL, UCL)
}

#--------------------------------------------------------------------------------
# Fit MR_weighted median on simulated data
#--------------------------------------------------------------------------------
fit_mr_median <- function(D, snps = 1:D$p){
                mr_weighted_median(D$beta_hat_x[snps], 
                                   D$beta_hat_y[snps], 
                                   sqrt(D$sigma2_x[snps]), 
                                   sqrt(D$sigma2_y[snps])) %>% 
                as_tibble() %>% 
                mutate(Estimate = b,
                       StdError = se,
                       LCL = Estimate - qnorm(0.975) * StdError,
                       UCL = Estimate + qnorm(0.975) * StdError
                ) %>%
                select(Estimate, StdError, LCL, UCL)
}

#--------------------------------------------------------------------------------
# Fit MR_Egger on simulated data
#--------------------------------------------------------------------------------
fit_mr_egger <- function(D, snps = 1:D$p){
                 mr_egger_regression(D$beta_hat_x[snps], 
                                     D$beta_hat_y[snps], 
                                     sqrt(D$sigma2_x[snps]), 
                                     sqrt(D$sigma2_y[snps])
                  )$mod$coefficients %>% 
                  as_tibble() %>%
                  filter(row_number() == 2) %>%
                  mutate(StdError = `Std. Error`,
                         LCL = Estimate - qnorm(0.975) * StdError,
                         UCL = Estimate + qnorm(0.975) * StdError
                  ) %>%
                  select(Estimate, StdError, LCL, UCL)
                  
}

#-------------------------------------------
# Fit Constrained ML with model averaging
#------------------------------------------
fit_mr_cML <- function(D){
  res_mr_cML_ <- MRcML::mr_cML(D$beta_hat_x, D$beta_hat_y, sqrt(D$sigma2_x), sqrt(D$sigma2_y), n = D$n)
  tibble(
    Estimate = res_mr_cML_$MA_BIC_theta,
    StdError = res_mr_cML_$MA_BIC_se
  ) %>% mutate(
    LCL = Estimate - qnorm(0.975) * StdError,
    UCL = Estimate + qnorm(0.975) * StdError
  )
}

#---------------------------------------------------------------
# Fit Constrained ML with model averaging and data perturbation
#---------------------------------------------------------------
fit_mr_cML_DP <- function(D){
  res_mr_cML_DP_ <- mr_cML_DP(D$beta_hat_x, D$beta_hat_y, sqrt(D$sigma2_x), sqrt(D$sigma2_y), n = D$n)
  tibble(
    Estimate = res_mr_cML_DP_$MA_BIC_DP_theta,
    StdError = res_mr_cML_DP_$MA_BIC_DP_se
  ) %>% mutate(
    LCL = Estimate - qnorm(0.975) * StdError,
    UCL = Estimate + qnorm(0.975) * StdError
  )
}

#---------------------------
# Fit Quantile Regression
#---------------------------
fit_mr_quantile <- function(D, snps = 1:D$p, se = "boot", verbose = FALSE, tau = NULL, nboot = 1000){
  mr_quantile(D$beta_hat_y[snps], D$beta_hat_x[snps], D$sigma2_y[snps], D$sigma2_x[snps], 
              se = se, verbose = verbose, tau = tau, nboot = nboot)$coefficients %>%
    as_tibble()
}

fit_mr_quantile_mv <- function(D, se = "boot", verbose = FALSE, tau = NULL, nboot = 1000){
  mr_quantile_mv(D$beta_hat_y, as.matrix(D$beta_hat_x), D$sigma2_y, as.matrix(D$sigma2_x), 
                 se = se, verbose = verbose, tau = tau, nboot = nboot)$coefficients %>%
    as_tibble()
}

fit_mr_quantile_ma <- function(D, se = "boot", weights = "bic",
                               tau_seq = seq(0.1, 0.9, by = 0.05), 
                               nboot = 200){
  
  res <- lapply(tau_seq, function(tau_val) {
    fit <- mr_quantile(D$beta_hat_y, 
                       D$beta_hat_x,
                       D$sigma2_y,
                       D$sigma2_x,
                       se = se, 
                       tau = tau_val
    )
    
    tibble(tau = tau_val, fit$coefficients, aic = fit$aic, bic = fit$bic)
  }) %>%
    bind_rows()
  
  if (weights == "bic"){
    w_k <- exp(-res$bic / 2) / sum(exp(-res$bic / 2))
  } else if (weights == "aic"){
    w_k <- exp(-res$aic / 2) / sum(exp(-res$aic / 2))
  }
  
  theta_k <- res$Estimate
  theta_hat <- sum(w_k * theta_k)
  theta_se <- sum(w_k * sqrt(res$StdError^2 + (theta_k - theta_hat)^2))
  tibble(
    Estimate = theta_hat,
    StdError = theta_se,
    LCL = theta_hat - 1.96 * theta_se,
    UCL = theta_hat + 1.96 * theta_se
  )
}

#---------------------------------------------
# Fit MR_weighted Mode on simulated data
#---------------------------------------------
fit_mr_mode <- function(D, snps = 1:D$p){
  mr_weighted_mode(D$beta_hat_x[snps], 
                   D$beta_hat_y[snps], 
                   sqrt(D$sigma2_x[snps]), 
                   sqrt(D$sigma2_y[snps])) %>% 
    as_tibble() %>% 
    mutate(Estimate = b,
           StdError = se,
           LCL = Estimate - qnorm(0.975) * StdError,
           UCL = Estimate + qnorm(0.975) * StdError
    ) %>%
    select(Estimate, StdError, LCL, UCL)
}

#---------------------------------------------
# Fit MR_RAPS on simulated data
#---------------------------------------------
fit_mr_raps <- function(D, snps = 1:D$p){
  mr_raps(D$beta_hat_x[snps], 
          D$beta_hat_y[snps], 
          sqrt(D$sigma2_x[snps]), 
          sqrt(D$sigma2_y[snps])
  ) %>% 
    as_tibble() %>% 
    mutate(Estimate = b,
           StdError = se,
           LCL = Estimate - qnorm(0.975) * StdError,
           UCL = Estimate + qnorm(0.975) * StdError
    ) %>%
    select(Estimate, StdError, LCL, UCL)
}

#---------------------------------------------
# Fit MR_Lasso on simulated data
#---------------------------------------------
fit_mr_lasso <- function(D) {
  tryCatch(
    {
      res_mr_lasso <- mr_lasso(mr_input(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), 
                                       by = D$beta_hat_y, byse = sqrt(D$sigma2_y)))
      tibble(
        Estimate = res_mr_lasso$Estimate,
        StdError = res_mr_lasso$StdError
      ) %>% mutate(
        LCL = res_mr_lasso$CILower,
        UCL = res_mr_lasso$CIUpper
      )
    },
    error = function(e) {
      warning(paste("Failed:", e$message))
      return(NULL)
    }
  )
}

#---------------------------------------------
# Fit MR_PRESSO on simulated data
#---------------------------------------------
fit_mr_presso <- function(D){
  res_mr_presso <- MRPRESSO::mr_presso(
                             BetaOutcome = "beta_hat_y", 
                             BetaExposure = "beta_hat_x", 
                             SdOutcome = "sqrt_sigma2_y", 
                             SdExposure = "sqrt_sigma2_x",
                             OUTLIERtest = TRUE, 
                             DISTORTIONtest = TRUE,
                             data = data.frame(beta_hat_x = D$beta_hat_x,
                                               beta_hat_y = D$beta_hat_y,
                                               sqrt_sigma2_x = sqrt(D$sigma2_x),
                                               sqrt_sigma2_y = sqrt(D$sigma2_y)
                             )
  )
  
  p_global <- res_mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue
  
  if(!is.na(p_global) && p_global < 0.05){
      tibble(
        Estimate = res_mr_presso$`Main MR results`$`Causal Estimate`[2],
        StdError = res_mr_presso$`Main MR results`$`Sd`[2]
      ) %>% mutate(
        LCL = Estimate - qnorm(0.975) * StdError,
        UCL = Estimate + qnorm(0.975) * StdError,
        Pvalue = p_global
      )
  } else{
      tibble(
        Estimate = res_mr_presso$`Main MR results`$`Causal Estimate`[1],
        StdError = res_mr_presso$`Main MR results`$`Sd`[1]
      ) %>% mutate(
        LCL = Estimate - qnorm(0.975) * StdError,
        UCL = Estimate + qnorm(0.975) * StdError,
        Pvalue = p_global
      )
  }
}

#---------------------------------------------
# Fit MR-ContMix on simulated data
#---------------------------------------------
fit_mr_contmix <- function(D){
  res_mr_contmix <- mr_conmix(mr_input(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), 
                                      by = D$beta_hat_y, byse = sqrt(D$sigma2_y)))
    
    # Check for valid confidence intervals
    valid_ci <- which(
        res_mr_contmix$Estimate >= res_mr_contmix$CILower & 
        res_mr_contmix$Estimate <= res_mr_contmix$CIUpper
    )
    
    if (length(res_mr_contmix$CILower) > 1 | length(res_mr_contmix$CIUpper) > 1) {
      warning("Multiple valid confidence intervals found in MR-ContMix. 
              Returning only one.")
    }
    
    # Output results
    tibble(
      Estimate = res_mr_contmix$Estimate,
      StdError = NA
    ) %>% mutate(
      LCL = res_mr_contmix$CILower[valid_ci],
      UCL = res_mr_contmix$CIUpper[valid_ci],
      Pvalue = res_mr_contmix$Pvalue
    )
}

#---------------------------------------------
# Fit MR-Mix on simulated data
#---------------------------------------------
fit_mr_mix <- function(D){
  res_mr_mix <- MRMix::MRMix(D$beta_hat_x, D$beta_hat_y, 
                             sqrt(D$sigma2_x), sqrt(D$sigma2_y)
  )
  
  tibble(
    Estimate = res_mr_mix$theta,
    StdError = res_mr_mix$SE_theta
  ) %>% mutate(
    LCL = Estimate - qnorm(0.975) * StdError,
    UCL = Estimate + qnorm(0.975) * StdError,
  )
}

#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
fit_mr_clust <- function(D){
    ratios <- as.vector(D$beta_hat_y / D$beta_hat_x)
    ratios.se <- sqrt(D$sigma2_y / D$beta_hat_x^2 + D$beta_hat_y^2 * D$sigma2_x  / D$beta_hat_x^4)
    # ratios.se <- as.vector(sqrt(D$sigma2_y) / abs(D$beta_hat_x))
    
    res_em <- mr_clust_em(theta = ratios, 
                          theta_se = ratios.se, 
                          bx = D$beta_hat_x, 
                          by = D$beta_hat_y, 
                          bxse = sqrt(D$sigma2_x), 
                          byse = sqrt(D$sigma2_y),
                          obs_names = paste0("SNP", 1:D$p)
                         )
    
    res_em_best <- res_em$results$best
    table(res_em_best$cluster_class, res_em_best$cluster_mean)
    
    # Find cluster with the most number of SNPs
    nclust <- length(unique(res_em_best$cluster))
    main_cluster <- names(which.max(table(res_em_best$cluster)))
    snps <- filter(res_em_best, cluster == as.numeric(main_cluster)) %>%
            mutate(rows = substr(observation, 4, 6)) %>%
            pull(rows) %>%
            as.numeric()
    
    # Fit IVW estimator on main cluster
    (if (nclust == 1){
      theta_main <- fit_mr_ivw(D, snps)
    } else {
      cluster_summary <- distinct(arrange(res_em_best, cluster), cluster, pi_hat, cluster_mean)
      theta_main <- data.frame(Estimate = cluster_summary$cluster_mean[as.numeric(main_cluster)],
                                StdError = sqrt(1 / colSums(res_em$r_ij * ratios.se^-2))[as.numeric(main_cluster)]
      ) %>%
          mutate(LCL = Estimate - qnorm(0.975) * StdError, UCL = Estimate + qnorm(0.975) * StdError) %>%
          as.matrix()
    })
    
    return(
      cbind(
           theta_main,
           fit_mr_median(D, snps),
           nclust = nclust,
           qhat = 1 - length(snps) / D$p,
           adj.rand = adjustedRandIndex(res_em_best$cluster_mean, D$mu),
           rand = rand.index(res_em_best$cluster_mean, D$mu)
           )
      )
}

#--------------------------------------------------------------------------------
# Fit Laplace mixture on simulated data
#--------------------------------------------------------------------------------
fit_mr_laplace <- function(D, dist = "laplace", k = NULL, tau = 0.5, null_cluster = FALSE, junk_cluster = FALSE){
  
  # Compute ratio estimates and their standard errors
  ratios <- D$beta_hat_y / D$beta_hat_x
  ratios.se <- sqrt(D$sigma2_y / D$beta_hat_x^2 + D$beta_hat_y^2 * D$sigma2_x  / D$beta_hat_x^4)
  # ratios.se <- sqrt(D$sigma2_y / D$beta_hat_x^2)
  
  # Fit model
  res_em <- laplace_mix_mr(dist, ratios, ratios.se, k, tau, null_cluster = null_cluster, junk_cluster = junk_cluster)
  table(res_em$results$cluster_class, res_em$results$cluster_median)
  
  # Find cluster with the most number of SNPs
  nclust <- length(unique(res_em$results$cluster))
  main_cluster <- names(which.max(table(res_em$results$cluster)))
  snps <- filter(res_em$results, cluster == main_cluster) %>%
    mutate(rows = substr(SNP, 4, 6)) %>%
    pull(rows) %>%
    as.numeric()
  
  # Fit weighted likelihood bootstrap
  Estimate <- first(filter(res_em$results, cluster == main_cluster))$cluster_median
  wlbs <- mr_laplace_boot(res_em, nboot = 200)
  boot_se <- unname(wlbs[as.numeric(main_cluster), "se_mu"])
    
  return(
    cbind(
      Estimate = Estimate,
      StdError = boot_se,
      LCL = Estimate - qnorm(0.975) * boot_se,
      UCL = Estimate + qnorm(0.975) * boot_se,
      fit_mr_median(D, snps),
      nclust = nclust,
      qhat = 1 - length(snps) / D$p,
      adj.rand = adjustedRandIndex(res_em$results$cluster_median, D$eta + D$theta),
      rand = rand.index(res_em$results$cluster_median, D$eta + D$theta)
    )
  )

}

# Check function
rho_tau <- function(u, tau) u * ( tau - (u < 0) )

# Weighted quantile function
weighted_quantile <- function(x, w, prob = 0.5) {
  
  # sort
  o <- order(x)
  x <- x[o]; w <- w[o]
  
  # cumulative weight
  cw <- cumsum(w) - 0.5 * w
  cw <- cw / sum(w)
  
  # find indices where cw < prob
  idx <- which(cw < prob)
  
  if (length(idx) == 0) return(first(x))         # prob is smaller than all cw
  if (length(idx) == length(x)) return(last(x)) # prob is larger than all cw
  
  k <- max(idx)
  
  # linear interpolation
  ratio <- (prob - cw[k]) / (cw[k + 1] - cw[k])
  x[k] + (x[k + 1] - x[k]) * ratio
}

# EM algorithm for mixture with heteroscedastic scales
laplace_mix_em <- function(dist, ratios, ratios.se, K, tau, weights, max_iter, tol, null_cluster, junk_cluster) {
  
  # ---- Initialization ----
  n <- length(ratios)
  
  if (dist == "laplace"){
    lambda <- 1
    w_i <- 1 / ratios.se
    lambda_i <- lambda * w_i
    isnull_tau <- is.null(tau)
    tau <- ifelse(isnull_tau, runif(1, 0.1, 0.9), tau)
  } else if (dist == "gaussian"){
    lambda <- NULL
    sigma_i <- ratios.se
  }
  
  # Initial cluster assignments
  if (K > 0){
    init_clust <- cluster::pam(ratios, K, diss = FALSE, metric = "manhattan")
    cluster <- init_clust$clustering
    mu <- c(init_clust$medoids); tmp_mu <- mu

    # Initial mixture proportions
    pi <- rep(1e-12, K)
    pi[as.numeric(names(table(cluster)))] <- table(cluster) / n
  } else {
    mu <- 0
    tmp_mu <- 0
    pi <- 1
  }
  
  # Null cluster
  if (K > 0 & null_cluster) {
    mu[K+1] <- 0
    pi_null <- runif(1, 0.05, 0.4)
    pi <- c( pi * ( 1 - pi_null ), pi_null )
  }
  
  # Junk cluster
  if (junk_cluster){
    mu[K+null_cluster+1] <- median(ratios)
    b_junk <- diff(range(ratios)) + 2 * max(ratios.se)
    pi_junk <- runif(1, 0.05, 0.4)
    pi <- c( pi * ( 1 - pi_junk ), pi_junk )
  }
  
  # ---- EM algorithm ----
  loglik_prev <- -Inf
  for (iter in 1:max_iter) {
    
    # ---- E-step ----
    log_r <- matrix(NA, n, K + null_cluster + junk_cluster)
    for (k in seq_len(K + null_cluster)) {
      if (dist == "laplace") {
        log_dens <- log(tau * ( 1 - tau )) + log(lambda_i) - lambda_i * rho_tau(ratios - mu[k], tau)
      } else if (dist == "gaussian") {
        log_dens <- dnorm(ratios, mean = mu[k], sd = sigma_i, log = TRUE)
      }
      log_r[, k] <- log(pi[k]) + log_dens
    }
    
    if (junk_cluster){
      log_junk_dens <- -log(2 * b_junk) - abs(ratios - mu[K + null_cluster + 1]) / b_junk
      log_r[, K + null_cluster + 1] <- log(pi[K + null_cluster + 1]) + log_junk_dens
    }
    
    # normalize rows
    log_row_sum <- apply(log_r, 1, logsumexp)
    r <- exp(log_r - log_row_sum)
    
    # ---- M-step ----
    # Update mixture proportions
    pi <- pmax(colSums(weights * r) / sum(weights), 1e-12)
    
    # Update tau
    if (isnull_tau & dist == "laplace"){
      a <- sum(weights * lambda_i * rowSums((r * outer(ratios, mu, "-"))[, seq_len(K + null_cluster), drop = FALSE]))
      c <- sum(weights * r[, seq_len(K + null_cluster)])
      roots <- ( a + 2 * c + c(-1, 1) * sqrt(a^2 + 4 * c^2)) / ( 2 * a )
      roots <- roots[is.finite(roots)]
      valid <- roots[roots > 0 & roots < 1]
      if (length(valid) >= 1) {
        # if two valid (unlikely), take the one closest to 0.5
        tau <- valid[which.min(abs(valid - 0.5))]
      } else {
        tau <- 0.5
      }
    }
    
    # Update lambda
    if (dist == "laplace"){
        lambda <- sum(weights * r[, seq_len(K + null_cluster)]) / sum( weights * w_i * rowSums((r * rho_tau(outer(ratios, mu, "-"), tau))[, seq_len(K + null_cluster), drop = FALSE]) )
        lambda <- min(lambda, sqrt(1 - 2 * tau * ( 1 - tau)) / (tau * (1 - tau)))
        lambda_i <- lambda * w_i
    }
    
    # Update cluster centers
    for (k in seq_len(K)) {
      if (dist == "laplace"){
          w <- weights * lambda_i * r[, k]
          tmp_mu[k] <- weighted_quantile(ratios, w, tau)
      } else if (dist == "gaussian"){
          w <- weights * r[, k] / (sigma_i^2)
          tmp_mu[k] <- sum(w * ratios) / sum(w)
      }
      mu[k] <- tmp_mu[k]
    }
    
    if (null_cluster){
      tmp_mu[K+1] <- weighted_quantile(ratios, weights * lambda_i * r[, K + 1], tau)
    }
    
    if (junk_cluster){ 
      tmp_mu[K + null_cluster + 1] <- weighted_quantile(ratios, weights * r[, K + null_cluster + 1] / b_junk, 0.5)
    }
    
    # ---- Log-likelihood ----
    dens_mat <- sapply(seq_len(K + null_cluster), function(k) {
      if (dist == "laplace") {
        pi[k] * lambda_i * tau * (1 - tau) * exp( -lambda_i * rho_tau(ratios - mu[k], tau) )
      } else if (dist == "gaussian") {
        pi[k] * dnorm(ratios, mean = mu[k], sd = sigma_i)
      }
    })
    
    if (junk_cluster) {
      junk_dens <- pi[K+null_cluster+1] / ( 2 * b_junk ) * exp( -abs( ratios - mu[K+null_cluster+1] ) / b_junk)
      dens_mat <- cbind(dens_mat, junk_dens)
    }
    
    loglik <- sum(weights * matrixStats::rowLogSumExps(log(dens_mat)))
    
    # Check convergence
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  if (iter == max_iter) {
    warning(paste("EM algorithm did not converge within the maximum number of iterations for K =", K))
  }
  
  # ---- Output ----
  cluster_class <- apply(r, 1, which.max)
  if (null_cluster) cluster_class[cluster_class == (K + 1)] <- "Null"
  if (junk_cluster) cluster_class[cluster_class == (K + null_cluster + 1)] <- "Junk"
    
  list(
    mu = tmp_mu,
    pi = pi,
    r = r,
    probability = apply(r, 1, max),
    loglik = loglik,
    iterations = iter,
    cluster = apply(r, 1, which.max),
    cluster_class = as.character(cluster_class),
    tau = tau,
    lambda = lambda
  )
}

# Fit EM for a range of K = 1 to p and select best K using BIC
laplace_mix_mr <- function(dist, ratios, ratios.se, k = NULL, tau = NULL, weights = rep(1, length(ratios)), max_iter = 5000, tol = 1e-5, 
                           snps = paste0("SNP", 1:length(ratios)), stop_bic_iter = 5, min_clust_search = 10, 
                           null_cluster = FALSE, junk_cluster = FALSE, nstarts = 5) {
  
  p <- length(ratios)
  K.max <- p - 1
  
  # If user provides a fixed number of clusters
  if (!is.null(k)) {
    fit <- laplace_mix_em(dist, ratios, ratios.se, k, tau, weights, max_iter, tol, null_cluster, junk_cluster)
    df <- ( length(fit$pi) - 1 ) + k # ( mixture proportions - 1 ) + cluster centers
    bic <- -2 * fit$loglik + df * log(p) # - 2 * sum(fit$r * log(fit$r))
    
    return(list(
      results = data.frame(
        SNP = snps,
        cluster = fit$cluster,
        pi_hat = round(fit$pi[fit$cluster], 3),
        probability = round(fit$probability, 3),
        cluster_median = round(fit$mu[fit$cluster], 3),
        theta = ratios,
        theta_se = ratios.se,
        cluster_class = fit$cluster_class
      ),
      pi = fit$pi,
      mu = fit$mu,
      r = fit$r,
      bic = bic,
      dist = dist,
      tau = fit$tau,
      lambda = fit$lambda,
      junk_cluster = junk_cluster
    ))
  }
  
  # Else search for best K using BIC
  best_bic <- Inf
  best_fit <- NULL
  bic <- numeric(K.max)  # store BIC for each K
    
  for (K in (1-null_cluster):K.max) {
    bic_temp <- numeric(nstarts)
    for (j in 1:nstarts){
        fit <- laplace_mix_em(dist, ratios, ratios.se, K, tau, weights, max_iter, tol, null_cluster, junk_cluster)
        
        # Calculate BIC
        df <- ( length(fit$pi) - 1 ) + K  # ( mixture proportions - 1 ) + cluster centers
        bic_temp[j] <- -2 * fit$loglik + df * log(p) # - 2 * sum(fit$r * log(fit$r))
        
        if (bic_temp[j] < best_bic) {
          best_bic <- bic_temp[j]
          best_fit <- fit
          best_fit$K <- K
          best_fit$bic <- bic_temp[j]
        }
    }
    bic[K+null_cluster] <- min(bic_temp)
    
    # ---- Early stopping rule ----
    if (K >= min_clust_search) {
      last_vals <- bic[(K + null_cluster - stop_bic_iter + 1):(K + null_cluster)]
      prev_vals <- bic[(K + null_cluster - stop_bic_iter):(K + null_cluster - 1)]
      
      # if BIC has increased monotonically in the last stop_bic_iter steps
      if (all(last_vals > prev_vals)) {
        message("Early stopping at K = ", K, 
                " (BIC worsened for ", stop_bic_iter, " consecutive steps).")
        break
      }
    }
  }
  
  return(list(
      results = data.frame(
        SNP = snps,
        cluster = best_fit$cluster,
        pi_hat = round(best_fit$pi[best_fit$cluster], 3),
        probability = round(best_fit$probability, 3),
        cluster_median = round(best_fit$mu[best_fit$cluster], 3),
        theta = ratios,
        theta_se = ratios.se,
        cluster_class = best_fit$cluster_class
      ),
      pi = best_fit$pi,
      mu = best_fit$mu,
      r = best_fit$r,
      bic = best_fit$bic,
      dist = dist,
      tau = best_fit$tau,
      lambda = best_fit$lambda,
      junk_cluster = junk_cluster,
      bic_list = bic[seq_len(K+null_cluster)]
      )
  )
}

mr_laplace_boot <- function(fit, nboot = 1000, alpha = 1, max_iter=5000, tol = 1e-5){
  
  # Initialize parameters from initial fit
  dist <- fit$dist
  ratios <- fit$results$theta
  ratios.se <- fit$results$theta_se
  r <- fit$r
  n <- nrow(r)
  K <- ncol(r)
  tau <- fit$tau
  lambda <- fit$lambda
  w_i <- 1 / ratios.se
  lambda_i <- lambda * w_i
  mu <- numeric(K)
  b_junk <- diff(range(ratios)) + 2 * max(ratios.se)
  junk_cluster <- fit$junk_cluster
  
  # Weighted likelihood bootstrap
  boot_samples <- replicate(nboot, {
    weights <- rgamma(n, shape = alpha); weights <- weights / mean(weights)
      
      # ---- EM algorithm ----
      loglik_prev <- -Inf
      for (iter in 1:max_iter) {
        
        # ---- M-step ----
        # Update mixture proportions
        pi <- pmax(colSums(weights * r) / sum(weights), 1e-12)
        
        # Update cluster centers
        for (k in 1:( K - junk_cluster )) {
          if (dist == "laplace"){
              w <- weights * lambda_i * r[, k]
              mu[k] <- weighted_quantile(ratios, w, tau)
          } else if (dist == "gaussian"){
              w <- weights * r[, k] / (sigma_i^2)
              mu[k] <- sum(w * ratios) / sum(w)
          }
        }
        
        # Update lambda
        if (dist == "laplace"){
          lambda <- sum(weights * r[, 1:K]) / sum( weights * w_i * rowSums((r * rho_tau(outer(ratios, mu, "-"), tau))[, 1:K, drop = FALSE]) )
          lambda <- min(lambda, sqrt(1 - 2 * tau * ( 1 - tau)) / (tau * (1 - tau)))
          lambda_i <- lambda * w_i
        }
        
        # Update location and scale parameters for the junk cluster
        if (junk_cluster){ 
          mu[K] <- weighted_quantile(ratios, weights * r[, K + 1] / b_junk, 0.5)
         # b_junk <- sum( r[, K] * abs( ratios - mu[K] ) ) / sum( r[,] )
        }
        
        # ---- E-step ----
        log_r <- matrix(NA, n, K)
        for (k in 1:( K - junk_cluster )) {
          if (dist == "laplace") {
            log_dens <- log(tau * ( 1 - tau )) + log(lambda_i) - lambda_i * rho_tau(ratios - mu[k], tau)
          } else if (dist == "gaussian") {
            log_dens <- dnorm(ratios, mean = mu[k], sd = sigma_i, log = TRUE)
          }
          log_r[, k] <- log(pi[k]) + log_dens
        }
        
        if (junk_cluster){
          log_junk_dens <- -log(2 * b_junk) - abs(ratios - mu[K]) / b_junk
          log_r[, K] <- log(pi[K]) + log_junk_dens
        }
        
        # normalize rows
        log_row_sum <- apply(log_r, 1, logsumexp)
        r <- exp(log_r - log_row_sum)
        
        # ---- Log-likelihood ----
        dens_mat <- sapply(1:( K - junk_cluster ), function(k) {
          if (dist == "laplace") {
            pi[k] * lambda_i * tau * (1 - tau) * exp( -lambda_i * rho_tau(ratios - mu[k], tau) )
          } else if (dist == "gaussian") {
            pi[k] * dnorm(ratios, mean = mu[k], sd = sigma_i)
          }
        })
        
        if (junk_cluster) {
          junk_dens <- pi[K] / ( 2 * b_junk ) * exp( -abs( ratios - mu[K] ) / b_junk)
          dens_mat <- cbind(dens_mat, junk_dens)
        }
        
        loglik <- sum(weights * matrixStats::rowLogSumExps(log(dens_mat)))
        
        # Check convergence
        if (abs(loglik - loglik_prev) < tol) break
        loglik_prev <- loglik
      }
      
      if (iter == max_iter) {
        warning(paste("EM algorithm did not converge within the maximum number of iterations for K =", K))
      }
      
      # ---- Return parameters ----
      list(mu = mu, pi = pi)
      
  }, simplify = FALSE)
  
  se_mu <- do.call(rbind, lapply(boot_samples, `[[`, "mu")) %>% apply(2, sd)
  se_pi <- do.call(rbind, lapply(boot_samples, `[[`, "pi")) %>% apply(2, sd)
  
  return(cbind("se_mu" = se_mu,
               "se_pi" = se_pi
               )
         )
  
}

#-----------------------------------------------
# Multivariable quantile regression
#-----------------------------------------------
mr_quantile_mv <- function(beta_hat_y,
                        beta_hat_x,
                        sigma2_y,
                        sigma2_x,
                        tau = NULL,
                        lambda = NULL,
                        theta = NULL,
                        weights = rep(1, length(beta_hat_y)),
                        tol = 1e-5,
                        max_iter = 500,
                        se = "boot",
                        nboot = 1000,
                        verbose = FALSE) {
  
  # Dimensions
  n <- length(beta_hat_y)
  K <- ncol(beta_hat_x)
  
  # Initialize parameters
  theta_init <- theta
  if (is.null(theta_init)) {
    res   <- optimize(
      Q_stat,
      interval   = c(-10, 10),
      beta_hat_x = beta_hat_x,
      beta_hat_y = beta_hat_y,
      sigma2_x   = sigma2_x,
      sigma2_y   = sigma2_y
    )
    theta <- res$minimum
  } else {
    theta <- theta_init
  }
  w_i <- 1 / sqrt(sigma2_y + sigma2_x * theta^2)
  
  # w_i <- w_i * n / sum(w_i)  # normalize weights
  lambda_init <- lambda
  lambda <- ifelse(is.null(lambda_init), 1, lambda_init)
  lambda_i <- lambda * w_i

  tau_init <- tau
  tau <- ifelse(is.null(tau_init), 0.5, tau_init)
  loglik_prev <- -Inf
  
  # Optimization loop
  for (iter in 1:max_iter) {
    
    # Update tau
    if (is.null(tau_init)){
      a <- sum(weights * lambda_i * (beta_hat_y - beta_hat_x %*% theta))
      tau <- 0.5 - a / ( 2 * ( 2 * sum(weights) + sqrt ( a^2 + 4 * sum(weights)^2 ) ) )
      tau <- min(max(tau, 0.1), 0.9)
    }
    
    # Update theta
    theta <- rq(beta_hat_y ~ beta_hat_x - 1, weights = weights * w_i, tau = tau)$coefficients
    
    # Update lambda
    if (is.null(lambda_init)){
      lambda <- sum(weights) / sum(weights * w_i * rho_tau(beta_hat_y - beta_hat_x %*% theta, tau))
      lambda_i <- lambda * w_i
    }
    
    # Log-likelihood
    loglik <- sum(weights) * (log(lambda) + log(tau * (1 - tau))) -
      lambda * sum(weights * w_i * rho_tau(beta_hat_y - beta_hat_x %*% theta, tau))
    
    # Convergence
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  if (iter == max_iter)
    warning("Maximum iterations reached without convergence.")
  
  # Bootstrap SE
  if (is.null(se)){
    se <- NA
    
  } else if (se == "boot") {
    se <- replicate(nboot, {
      beta_hat_y_boot <- rnorm(n, mean = beta_hat_y, sd = sqrt(sigma2_y))
      beta_hat_x_boot <- sapply(1:K, 
                                function(k){
                                  rnorm(n, mean = beta_hat_x[, k], sd = sqrt(sigma2_x[, k]))
                                }
                        )
      mr_quantile_mv(
        beta_hat_y_boot,
        beta_hat_x_boot,
        sigma2_y,
        sigma2_x, 
        tau = tau_init,
        lambda = lambda_init,
        theta = theta,
        se = NULL
      )$coefficients$Estimate
    }) %>% sd()
  }
  
  # Print results if verbose
  if (verbose) {
    cat("Converged after", iter, "iterations\n")
    cat("Estimated tau:", tau, "\n")
    cat("Estimated lambda:", lambda, "\n")
  }
  
  # Return data frame
  list(
    coefficients = data.frame(
      Estimate = theta,
      StdError = se,
      LCL = theta - qnorm(0.975) * se,
      UCL = theta + qnorm(0.975) * se
    ),
    tau = tau,
    lambda = lambda,
    aic = -2 * loglik + 2 * (2 + !is.null(tau_init)),
    bic = -2 * loglik + log(n) * (2 + !is.null(tau_init))
  )
}

Q_stat <- function(theta, beta_hat_x, beta_hat_y, sigma2_x, sigma2_y) {
  w_j <- beta_hat_x^2 / (theta^2 * sigma2_x + sigma2_y)
  sum(w_j * (beta_hat_y - theta * beta_hat_x)^2)
}
