#--------------------------------------------------------------------------------
# Simulate summary-level data
#--------------------------------------------------------------------------------

# --- Laplace sampler ---
rLaplace <- function(n, mu = 0, b = 1) {
  u <- runif(n, -0.5, 0.5)
  mu - b * sign(u) * log(1 - 2 * abs(u))
}

sim_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation parameters
  #--------------------------------------------------------------------------------
  # args <- list(theta = 0.2, p = 100000, h2_x = 0.3, h2_y = 0.3, n = 50000, q = 0.6, m = 100, pi = 0.2)
  n <- args$n
  theta <- args$theta
  p <- args$p
  h2_x <- args$h2_x
  h2_y <- args$h2_y
  q <- args$q
  m <- args$m
  pi <- args$pi
  
  # Simulate p SNPs with m non-null effects, m non-null uncorrelated pleiotropic effects, 
  # pi is the proportion of the SNPs with both beta_x and alpha different from zero
  maf <- runif(p, 0.05, 0.5)
  
  id.x <- sample(1:p, m) %>% sort()
  id.alpha <- c(sample(id.x, pi * m), sample(setdiff(1:p, id.x), (1 - pi) * m)) %>% sort()
  id.all <- sort(unique(c(id.x, id.alpha)))
  
  sd.x <- 1 / sqrt(2 * maf[id.x] * (1 - maf[id.x]))
  sd.alpha <- 1 / sqrt(2 * maf[id.alpha] * (1 - maf[id.alpha]))
  
  # Simulate direct effects
  sigma2_x <- h2_x / m
  beta_x <- rep(0, p)
  beta_x[id.x] <- rnorm(m, mean = 0, sd = sqrt(sigma2_x)) * sd.x
  sum(2 * maf * (1 - maf) * beta_x^2)
  
  # Simulate uncorrelated pleiotropic effects and m * q correlated pleiotropic effects
  # if (q == 0 ){
  #     etavals <- rep(0, 4)
  # } else if (q == 0.2){
  #     etavals <- 0.85 + 0.015 * (theta == 0)
  # } else if (q == 0.4){
  #     etavals <- 0.6 + 0.0125 * (theta == 0)
  # } else if (q == 0.6){
  #     etavals <- 0.532 + 0.008 * (theta == 0)
  # }
  
  if (q == 0 ){
    etavals <- rep(0, 4)
  } else {
    etavals <- 0.6
  }
  
  etavals <- etavals * c(0, -1, 1, 2)
  eta <- rep(0, p)
  eta[id.all] <- rep(etavals, times = c( round((1 - q) * length(id.all)), rep( q * length(id.all) / 3, 3 ) ) )
  etasq <- etavals^2 %*% c(1 - q, rep(q / 3, 3))
  
  sigma2_alpha <- (h2_y - (theta^2 + etasq) * h2_x) / m
  alpha <- rep(0, p)
  # alpha[id.alpha] <- rnorm(m, mean = 0, sd = sqrt(sigma2_alpha)) * sd.alpha
  alpha[id.alpha] <- rLaplace(m, mu = 0, b = sqrt(rep(sigma2_alpha, m) / 2)) * sd.alpha
  sum(2 * maf * (1 - maf) * alpha^2)
  
  # Simulate outcome
  beta_y <- theta * beta_x + eta * beta_x + alpha
  sum(2 * maf * (1 - maf) * beta_y^2)
  
  # Generate summary statistics
  sigma2 <-  1 / ( 2 * n * maf * ( 1 - maf ) ) 
  beta_hat_x <- rnorm( p, beta_x, sqrt( sigma2 ) )
  beta_hat_y <- rnorm( p, beta_y, sqrt( sigma2 ) )
  
  # Select SNPs with p-value < 5e-8
  pval_x <- 2 * pnorm(-abs(beta_hat_x / sqrt(sigma2)))
  snps <- which(pval_x < 5e-8)
  
  # Compute heritabilities
  (tau2 <- theta^2 * h2_x / h2_y)
  (omega2 <- etasq * h2_x / h2_y)
    
  return(list("beta_hat_x" = beta_hat_x[snps], 
              "sigma2_x" = sigma2[snps], 
              "beta_hat_y" = beta_hat_y[snps], 
              "sigma2_y" = sigma2[snps], 
              "theta" = theta,
              "alpha" = alpha[snps],
              "eta" = eta[snps],
              "maf" = maf[snps],
              "n" = n, 
              "q" = q,
              "p" = length(snps)
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
  # args <- list(theta_1 = 0, p = 90, n = 10000, q = 0.6, m = 10, k = 4, rho = 0.25)
  n <- args$n 
  theta_1 <- args$theta_1
  p <- args$p
  q <- args$q
  m <- args$m
  k <- args$k
  rho <- args$rho
  
  # Simulate genetic data and confounders
  maf <- runif(p, 0.1, 0.3)
  G <- sapply(maf, function(maf) rbinom(2 * n, 2, maf))
  epsilon_u <- mvrnorm(2 * n, mu = rep(0, k), Sigma = diag(1, k))
  epsilon_x <- mvrnorm(2 * n, mu = rep(0, m), Sigma = diag(1 - rho, m) + rho)
  epsilon_y <- rnorm(2 * n)
  
  # Simulate latent factors
  cluster <- rep(1:k, times = c(p - q * p, rep(q * p/(k - 1), k - 1)))
  phi <- runif(p, -0.1, 0.1) * sapply(1:k, function(i) cluster == i)
  phi[, 1] <- 0
  U <- G %*% phi + epsilon_u
  # apply(G %*% phi, 2, var) / apply(U, 2, var)
  
  # Simulate exposure
  gamma <- replicate(m, runif(p, 0.1, 0.2) * sample(c(-1, 1), size = p, replace = TRUE))
  gamma[cluster != 1, ] <- 0
  delta_x <- matrix(runif(k*m, max = 1, min = -1), nrow = k)
  X <- G %*% gamma + U %*% delta_x + epsilon_x
  # apply(G %*% gamma, 2, var) / apply(X, 2, var)
  
  # Simulate outcome
  theta <- c(theta_1, -0.1, 0.1, rep(0, m - 3))
  alpha <- rep(0, p)
  alpha[sample(1:p, q * p)] <- runif(q * p, 0.2, 0.3)
  delta_y <- runif(k, max = 1, min = -1)
  Y <- X %*% theta + G %*% alpha + U %*% delta_y + epsilon_y
  
  # Summary stats for X
  sampX <- seq(1, n)
  regX <- apply(X[sampX, ], 2, function(xcol) sumstat(xcol, G[sampX, ]))
  beta_hat_x <- sapply(regX, `[[`, "beta_hat")
  sigma2_x    <- sapply(regX, `[[`, "sigma2")
    
  # Summary stats for Y
  sampY <- seq(n+1, 2*n)
  regY <- sumstat(Y[sampY], G[sampY, ])
  beta_hat_y <- regY$beta_hat
  sigma2_y <- regY$sigma2
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2_x, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2_y,
              "cluster" = cluster,
              "theta" = theta,
              "alpha" = alpha,
              "eta" = c(0, delta_y[-1]),
              "n" = n,
              "p" = p,
              "q" = q,
              "maf" = maf
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
# theta_ivw <- function(D, snps = 1:D$p){
# 
#   total.effect <- lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1)$coef[1]
#   resid.std.error <-summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$sigma
#   se.total.effect <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$coef[1,2]
#   ci.upper.total <- max(total.effect + qnorm(0.975) * se.total.effect / resid.std.error,
#                        total.effect + qt(0.975, df=length(snps)-1) * se.total.effect)
#   ci.lower.total <- min(total.effect - qnorm(0.975) * se.total.effect / resid.std.error,
#                        total.effect - qt(0.975, df=length(snps)-1) * se.total.effect)
# 
#   return(cbind("Estimate" = total.effect, "StdError" = se.total.effect, "LCL" = ci.lower.total, "UCL" = ci.upper.total))
# }
# 
# theta_mvivw <- function(D, snps = 1:D$p){
# 
#   direct.effect <- lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] + D$beta_hat_m[snps] - 1, weights = D$sigma2_y[snps]^-1)$coef[1]
#   resid.std.error <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] - 1, weights = D$sigma2_y[snps]^-1))$sigma
#   se.direct.effect <- summary(lm(D$beta_hat_y[snps] ~ D$beta_hat_x[snps] + D$beta_hat_m[snps] - 1, weights = D$sigma2_y[snps]^-1))$coef[1,2]
#   ci.upper.direct <- max(direct.effect + qnorm(0.975) * se.direct.effect / resid.std.error,
#                        direct.effect + qt(0.975, df=length(snps)-1) * se.direct.effect)
#   ci.lower.direct <- min(direct.effect - qnorm(0.975) * se.direct.effect / resid.std.error,
#                        direct.effect - qt(0.975, df=length(snps)-1) * se.direct.effect)
# 
#   return(cbind("Estimate" = direct.effect, "StdError" = se.direct.effect, "LCL" = ci.lower.direct, "UCL" = ci.upper.direct))
# }
# 
# fit_ivw <- function(D, snps = 1:D$p){
#   theta_T <- theta_ivw(D, snps)
#   theta_D <- theta_mvivw(D, snps)
#   theta <- cbind(
#           theta_T_hat = theta_T[,"Estimate"], theta_T_se = theta_T[,"StdError"],
#           theta_T_lcl = theta_T[,"LCL"], theta_T_ucl = theta_T[,"UCL"],
#           theta_D_hat = theta_D[,"Estimate"], theta_D_se = theta_D[,"StdError"],
#           theta_D_lcl = theta_D[,"LCL"], theta_D_ucl = theta_D[,"UCL"]
#         )
#   row.names(theta) <- NULL
#   return(theta)
# }

theta_ivw <- function(D, snps = 1:D$p){
  MRInputObject = mr_input(bx = D$beta_hat_x[snps], bxse = sqrt(D$sigma2_x[snps]), by = D$beta_hat_y[snps], byse = sqrt(D$sigma2_y[snps]))
  MR_IVW <- mr_ivw(MRInputObject)
  return(cbind("Estimate" = MR_IVW$Estimate, 
               "StdError" = MR_IVW$StdError,
               "LCL" = MR_IVW$CILower,
               "UCL" = MR_IVW$CIUpper
               )
         )
}

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
    return(cbind("Estimate" = MR_Egger$Estimate, 
                 "StdError" = MR_Egger$StdError.Est, 
                 "LCL" = MR_Egger$CILower.Est,
                 "UCL" = MR_Egger$CIUpper.Est
                 )
           )
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
theta_median <- function(D, snps = 1:D$p, w = rep(1, length(snps))){
    MR_Median <- summary(rq(D$beta_hat_y[snps]~-1+D$beta_hat_x[snps], weights = w * D$sigma2_y[snps]^-1), se = "boot")$coef
    return(cbind("Estimate" = MR_Median[1], 
                 "StdError" = MR_Median[2], 
                 "LCL" = MR_Median[1] - qt(0.975, df = length(snps) - 1) * MR_Median[2],
                 "UCL" = MR_Median[1] + qt(0.975, df = length(snps) - 1) * MR_Median[2]
                 )
           )
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
  row.names(theta) <- NULL
  return(theta)
}

#--------------------------------------------------------------------------------
# Fit MR_Clust on simulated data
#--------------------------------------------------------------------------------
fit_mrclust <- function(D){
    ratios <- as.vector(D$beta_hat_y / D$beta_hat_x)
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
    nclust <- length(unique(res_em_best$cluster))
    main_cluster <- names(which.max(table(res_em_best$cluster)))
    snps <- filter(res_em_best, cluster == as.numeric(main_cluster)) %>%
            mutate(rows = substr(observation, 4, 6)) %>%
            pull(rows) %>%
            as.numeric()
    
    # Fit IVW estimator on main cluster
    (if (nclust == 1){
      theta_main <- theta_ivw(D, snps)
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
           theta_median(D, snps),
           nclust = nclust,
           qhat = 1 - length(snps) / D$p,
           adj.rand = adjustedRandIndex(res_em_best$cluster_mean, D$eta + D$theta),
           rand = rand.index(res_em_best$cluster_mean, D$eta + D$theta)
           )
      )
}

fit_mvmrclust <- function(D){
  fit <- lm(D$beta_hat_y ~ D$beta_hat_m, weights = D$sigma2_y^-1)
  resid.y <- residuals(fit)
  resid.se <- sqrt(D$sigma2_y * (1 - hatvalues(fit)))
  ratios <- as.vector(resid.y / D$beta_hat_x)
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
  
  # Fit median estimator on main cluster
  nclust <- length(unique(res_em_best$cluster))
  theta_main <- fit_median(D, snps)
  
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
# Fit Bayesian Laplace mixture on simulated data
#--------------------------------------------------------------------------------
fit_mrlaplace_bayes <- function(D){
  
  # Fit model
  res_em <- mr_laplace(D$beta_hat_y, D$beta_hat_x, sigma2_y = D$sigma2_y)
  table(res_em$results$cluster_median)
  
  # Find cluster with the most number of SNPs
  nclust <- length(unique(res_em$results$cluster))
  main_cluster <- names(which.max(table(res_em$results$cluster)))
  snps <- filter(res_em$results, cluster == main_cluster) %>%
    mutate(rows = substr(SNP, 4, 6)) %>%
    pull(rows) %>%
    as.numeric()
  
  return(
    cbind(
      Estimate = first(filter(res_em$results, cluster == main_cluster))$cluster_median,
      StdError = first(filter(res_em$results, cluster == main_cluster))$cluster_se,
      LCL = first(filter(res_em$results, cluster == main_cluster))$cluster_lcl,
      UCL = first(filter(res_em$results, cluster == main_cluster))$cluster_ucl,
      nclust = nclust,
      qhat = 1 - length(snps) / D$p,
      adj.rand = adjustedRandIndex(res_em$results$cluster_median, D$eta + D$theta),
      rand = rand.index(res_em$results$cluster_median, D$eta + D$theta)
    )
  )
  
}

#--------------------------------------------------------------------------------
# Fit Laplace mixture on simulated data
#--------------------------------------------------------------------------------
fit_mrlaplace <- function(D, tau = 0.5, null_cluster = FALSE){
  
  # Compute ratio estimates and their standard errors
  ratios <- D$beta_hat_y / D$beta_hat_x
  ratios.se <- sqrt(D$sigma2_y / D$beta_hat_x^2)
  
  # Fit model
  res_em <- laplace_mix_mr(ratios, ratios.se, tau = tau, null_cluster = null_cluster)
  table(res_em$results$cluster_median)
  
  # Find cluster with the most number of SNPs
  nclust <- length(unique(res_em$results$cluster))
  main_cluster <- names(which.max(table(res_em$results$cluster)))
  snps <- filter(res_em$results, cluster == main_cluster) %>%
    mutate(rows = substr(SNP, 4, 6)) %>%
    pull(rows) %>%
    as.numeric()
  
  return(
    cbind(
      theta_median(D, snps),
      nclust = nclust,
      qhat = 1 - length(snps) / D$p,
      adj.rand = adjustedRandIndex(res_em$results$cluster_median, D$eta + D$theta),
      rand = rand.index(res_em$results$cluster_median, D$eta + D$theta)
    )
  )

}

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

# Update tau using quadratic formula
update_tau <- function(ratios, r, mu, b_i, eps = 1e-12) {
  
  # Compute weighted average of distances for each cluster
  n <- length(ratios)
  K <- ncol(r)
  c <- sum(r * sapply(1:K, function(k) (ratios - mu[k]) / b_i)) / n
  if (abs(c) < eps) return(0.5)
  
  # Find the roots of the quadratic equation
  roots <- ( c + 2 + c(-1, 1) * sqrt(c^2 + 4)) / (2 * c)
  roots <- roots[is.finite(roots)]
  valid <- roots[roots > 0 & roots < 1]
  
  if (length(valid) >= 1) {
    # if two valid (unlikely), take the one closest to 0.5
    return(valid[which.min(abs(valid - 0.5))])
  } else {
    # numeric fallback
    return(0.5)
  }
}

# EM for Laplace mixture with heteroscedastic scales
laplace_mix_em <- function(ratios, ratios.se, tau, K, max_iter, tol, null_cluster) {
  
  # ---- Initialization ----
  n <- length(ratios)
  if (is.null(tau)) {
    tau <- 0.5
    is.null.tau <- TRUE
  } else {
    is.null.tau <- FALSE
  }
  b_i <- ratios.se * sqrt( tau * ( 1 - tau ) / ( 1 - 2 * tau + 2 * tau^2 ) )
  
  # Initial cluster assignments
  cluster <- cluster::pam(ratios, k = K, diss = FALSE, metric = "manhattan", cluster.only=TRUE)
  mu <- numeric(K); tmp_mu <- numeric(K)
  for (k in 1:K) {
    xk <- ratios[cluster == k]
    mu[k] <- weighted_quantile(xk, rep(1, length(xk)), prob = 0.5)
  }
  
  # Initial mixture proportions
  pi <- table(cluster) / n
  
  # check-loss function
  rho_tau <- function(u) u * (tau - (u < 0))
  
  # Determine which cluster should be null
  if (null_cluster) {
    null_idx <- which.min(abs(mu))  # Closest to zero
    mu[null_idx] <- 0
  }
  
  # ---- EM algorithm ----
  loglik_prev <- -Inf
  for (iter in 1:max_iter) {
    
    # ---- E-step ----
    r <- matrix(0, n, K)
    for (k in 1:K) {
      u <- (ratios - mu[k]) / b_i
      dens <- tau * (1 - tau) / b_i * exp(-rho_tau(u))
      r[, k] <- pi[k] * dens
    }
    
    # normalize rows
    r <- r / rowSums(r)
    
    # ---- M-step ----
    # Update mixture proportions
    pi <- colMeans(r)
    
    # Update tau and b_i's if tau is not fixed
    if (is.null.tau){
      tau <- update_tau(ratios, r, mu, b_i)
      b_i <- ratios.se * sqrt( tau * ( 1 - tau ) / ( 1 - 2 * tau + 2 * tau^2 ) )
    }
    
    # Update cluster centers
    for (k in 1:K) {
      w <- r[, k] / b_i
      tmp_mu[k] <- weighted_quantile(ratios, w, tau)
      if (null_cluster && k == null_idx) {
        next  # Skip updating the null cluster
      }
      mu[k] <- tmp_mu[k]
    }
    
    # ---- Log-likelihood ----
    loglik <- sum(log(rowSums(sapply(1:K, function(k) {
      u <- (ratios - mu[k]) / b_i
      pi[k] * tau * (1 - tau) / b_i * exp(-rho_tau(u))
    }))))
    
    # Check convergence
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }
  
  if (iter == max_iter) {
    warning(paste("EM algorithm did not converge within the maximum number of iterations for K =", K))
  }
  
  
  # ---- Output ----
  list(
    mu = tmp_mu,
    tau = tau,
    pi = pi,
    r = r,
    b_i = b_i,
    probability = apply(r, 1, max),
    loglik = loglik,
    iterations = iter,
    cluster = apply(r, 1, which.max)
  )
}

# Fit EM for a range of K = 1 to p and select best K using BIC
laplace_mix_mr <- function(ratios, ratios.se, k = NULL, tau = NULL, n = NULL, gamma = 0, max_iter = 5000, tol = 1e-5, 
                           snps = paste0("SNP", 1:length(ratios)), stop_bic_iter = 5, min_clust_search = 10, null_cluster = FALSE) {
  
  p <- length(ratios)
  K.max <- p - 1
  if (is.null(n)) { n <- p }
  
  # If user provides a fixed number of clusters
  if (!is.null(k)) {
    fit <- laplace_mix_em(ratios, ratios.se, tau, k, max_iter, tol, null_cluster)
    df <- 2 * k - 1
    bic <- -2 * fit$loglik + df * log(n) + 2 * gamma * p * log(k)
    
    return(list(
      results = data.frame(
        SNP = snps,
        cluster = fit$cluster,
        pi_hat = fit$pi,
        probability = fit$probability,
        cluster_quantile = fit$mu[fit$cluster],
        theta = ratios,
        theta_se = ratios.se
      ),
      tau = fit$tau,
      pi = fit$pi,
      mu = fit$mu,
      r = fit$r,
      b_i = fit$b_i,
      K = k,
      bic = bic
    ))
  }
  
  # Else search for best K using BIC
  best_bic <- Inf
  best_fit <- NULL
  bic <- numeric(K.max)  # store BIC for each K
  
  for (K in 1:K.max) {
    fit <- laplace_mix_em(ratios, ratios.se, tau, K, max_iter, tol, null_cluster)
    
    # Calculate BIC
    df <- 2 * K - 1  # mixture proportions + cluster centers
    bic[K] <- -2 * fit$loglik + df * log(n) + 2 * gamma * p * log(K)
    
    if (bic[K] < best_bic) {
      best_bic <- bic[K]
      best_fit <- fit
      best_fit$K <- K
      best_fit$bic <- bic[K]
    }
    
    # ---- Early stopping rule ----
    if (K >= min_clust_search) {
      last_vals <- bic[(K - stop_bic_iter + 1):K]
      prev_vals <- bic[(K - stop_bic_iter):(K - 1)]
      
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
      theta_se = ratios.se
    ),
    tau = best_fit$tau,
    pi = best_fit$pi,
    mu = best_fit$mu,
    r = best_fit$r,
    b_i = best_fit$b_i
  ))
}