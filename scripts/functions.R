#--------------------------------------------------------------------------------
# Simulate data
#--------------------------------------------------------------------------------
sim_dat <- function(args){
  
  #--------------------------------------------------------------------------------
  # Define simulation scenario, fraction of invalid IVs and theta
  #--------------------------------------------------------------------------------
  # args <- data.frame(scenario = 1, frac = 0.1, theta = 0)
  theta <- args$theta
  n <- 10000
  p <- 30
  
  #--------------------------------------------------------------------------------
  # Scenario 1 : Balanced pleiotropy, InSIDE assumption met
  #--------------------------------------------------------------------------------
  if (args$scenario == 1){
    phi <- rep(0, p) 
    alpha <- rnorm(p, 0, 0.15)
  }
  
  #--------------------------------------------------------------------------------
  # Scenario 2 : Directional pleiotropy, InSIDE assumption met
  #--------------------------------------------------------------------------------
  if (args$scenario == 2){
    phi <- rep(0, p) 
    alpha <- rnorm(p, 0.1, 0.075)
  }
  
  #--------------------------------------------------------------------------------
  # Scenario 3 : Positive directional pleiotropy, InSIDE assumption violated
  #--------------------------------------------------------------------------------
  if (args$scenario == 3){
    phi <- runif(p, 0, 0.1)
    alpha <- rnorm(p, 0.1, 0.075)
  }
  
  #--------------------------------------------------------------------------------
  # Scenario 4 : Negative directional pleiotropy, InSIDE assumption violated
  #--------------------------------------------------------------------------------
  if (args$scenario == 4){
    phi <- runif(p, 0, 0.1)
    alpha <- rnorm(p, -0.1, 0.075)
  }
  
  #--------------------------------------------------------------------------------
  # For each scenario either 30%, 50% or 70% of genetic variants are invalid
  #--------------------------------------------------------------------------------
  valid.IVs <- sample(1:p, (1-args$frac)*p)
  invalid.IVs <- setdiff(1:p, valid.IVs)
  phi[valid.IVs] <- 0
  alpha[valid.IVs] <- 0
  
  #--------------------------------------------------------------------------------
  # Simulate data
  #--------------------------------------------------------------------------------
  # Simulate genetic data and confounder
  maf <- runif(p, 0.1, 0.5)
  G <- sapply(1:p, function(j) rbinom(2 * n, 2, maf[j]))
  U <- G %*% phi + rnorm(2 * n)
  
  # Simulate uniform genetic effects such that the total proportion of variance explained 
  # in the exposure by direct effects of the genetic variants is on average 10%
  sigma2_gamma <- 1/9 * (var(U) + 1) * n / sum(G^2)
  gamma <- runif(p, 0, sqrt(12 * sigma2_gamma))
  
  # Simulate exposure and outcome
  X <- G %*% gamma + U + rnorm(2 * n)
  Y <- G %*% alpha + theta * X + U + rnorm(2 * n)
  
  #--------------------------------------------------------------------------------
  # Genetic variant-trait association
  #--------------------------------------------------------------------------------
  sampX = seq(1, n)
  sampY = seq((n+1), (2*n))
  
  # Summary stats for X
  beta_hat_x <- unname(sapply(1:p, function(i) coef(lm(X[sampX] ~ G[sampX, i]))[-1]))
  sigma2_x <- unname(sapply(1:p, function(i) coefficients(summary(lm(X[sampX] ~ G[sampX, i])))[2, "Std. Error"]^2))
  
  # Summary stats for Y
  beta_hat_y <- unname(sapply(1:p, function(i) coef(lm(Y[sampY] ~ G[sampY, i]))[-1]))
  sigma2_y <- unname(sapply(1:p, function(i) coefficients(summary(lm(Y[sampY] ~ G[sampY, i])))[2, "Std. Error"]^2))
  
  return(list("beta_hat_x" = beta_hat_x, 
              "sigma2_x" = sigma2_x, 
              "beta_hat_y" = beta_hat_y, 
              "sigma2_y" = sigma2_y, 
              "valid.IVs" = valid.IVs,
              "scenario" = args$scenario,
              "frac" = args$frac,
              "theta" = args$theta
              )
         )
}

#--------------------------------------------------------------------------------
# Fit MR_IVW on simulated data
#--------------------------------------------------------------------------------
theta_ivw <- function(D){
    MRInputObject = mr_mvinput(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), by = D$beta_hat_y, byse = sqrt(D$sigma2_y))
    MR_IVW <- mr_mvivw(MRInputObject)
    return(cbind("Estimate" = unname(MR_IVW$Estimate[1]), "StdError" = unname(MR_IVW$StdError[1]), "Pvalue" = unname(MR_IVW$Pvalue[1])))
}

#--------------------------------------------------------------------------------
# Fit MR_Egger on simulated data
#--------------------------------------------------------------------------------
theta_egger <- function(D){
    MRInputObject = mr_mvinput(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), by = D$beta_hat_y, byse = sqrt(D$sigma2_y))
    MR_Egger <- mr_mvegger(MRInputObject)
    return(cbind("Estimate" = unname(MR_Egger$Estimate[1]), "StdError" = unname(MR_Egger$StdError.Est[1]), "Pvalue" = unname(MR_Egger$Pvalue.Est[1])))
}

#--------------------------------------------------------------------------------
# Fit MR_Median on simulated data
#--------------------------------------------------------------------------------
theta_median <- function(D){
    p <- length(D$beta_hat_y)
    k <- dim(D$beta_hat_x)[2]
    theta_Median <- rq(D$beta_hat_y~-1+D$beta_hat_x, weights = D$sigma2_y^-1)$coefficients[1]
    names(theta_Median)[1] <- "Estimate"
    theta_Median["Std. Error"] <- replicate(1000, {
      beta_hat_y_k <- rnorm(p, D$beta_hat_y, sqrt(D$sigma2_y))
      beta_hat_x_k <- t(sapply(1:p, function(j) mvrnorm(1, D$beta_hat_x[j, ], diag(D$sigma2_x[j, ]))))
      rq(beta_hat_y_k~-1+beta_hat_x_k, weights = D$sigma2_y^-1)$coefficients[1]
    }) %>% sd()
    theta_Median["Pr(>|t|)"] <- 2*pt(q=theta_Median["Estimate"] / theta_Median["Std. Error"], df = p-k, lower.tail = FALSE)
    return(theta_Median)
}

#--------------------------------------------------------------------------------
# Fit MR_lasso on simulated data
#--------------------------------------------------------------------------------
theta_lasso <- function(D){
    p <- length(D$beta_hat_y)
    k <- dim(D$beta_hat_x)[2]
    
    # Change allele orientation and fit lasso model
    BYG <- D$beta_hat_y*sign(D$beta_hat_x)
    BXG <- abs(D$beta_hat_x)
    fit.lasso <- glmnet(cbind(diag(rep(1, p)), BXG), 
                        BYG, 
                        weights = D$sigma2_y^-1, 
                        intercept = FALSE,
                        penalty.factor = c(rep(1, p), rep(0, k))
    )
    
    # Select best model using the heterogeneity stopping rule
    nlambda <- fit.lasso$dim[2]
    lasso.path.IVs <- sapply(1:nlambda, function(i) which(fit.lasso$beta[1:p, i] == 0)) %>% rev()
    df.i <- lapply(1:nlambda, function(i) data.frame(D$beta_hat_y, D$beta_hat_x, D$sigma2_y^-1)[lasso.path.IVs[[i]], ])
    lm.i <- lapply(1:nlambda, function(i) summary(lm(D.beta_hat_y~-1+X1+X2+X3+X4, weights = D.sigma2_y..1, data = df.i[[i]])))
    J <- sapply(lasso.path.IVs, length)
    RSE <- sapply(1:nlambda, function(i) lm.i[[i]]$sigma)
    
    het <- which((RSE[-1] > 1) & (diff(RSE) > (qchisq(0.95, df = 1) / J[-1])))
    if (length(het) == 0){
      lam_pos = 1
    } else {
      lam_pos = min(het)
    }
    
    # Fit IVW estimator using the selected IVs
    v <- lasso.path.IVs[[lam_pos]]
    theta_lasso_ivw <- mr_mvivw(mr_mvinput(bx = D$beta_hat_x[v, ], 
                                           bxse = sqrt(D$sigma2_x)[v, ],
                                           by = D$beta_hat_y[v], 
                                           byse = sqrt(D$sigma2_y)[v]
    )
    )
    return(cbind("Estimate" = unname(theta_lasso_ivw$Estimate[1]), "StdError" = unname(theta_lasso_ivw$StdError[1]), "Pvalue" = unname(theta_lasso_ivw$Pvalue[1])))
}