library(quantreg)
library(MendelianRandomization)

# Load simulated data
load("../Data/sim_dat_null.R")
load("../Data/sim_dat.R")

# Load example
D <- D_S3_50[[1]]
n <- length(D$beta_hat_y)
D$beta_hat_y_Egger <- D$beta_hat_y * sign(D$beta_hat_x)
D$beta_hat_x_Egger <- D$beta_hat_x * sign(D$beta_hat_x)
MRInputObject = mr_input(bx = D$beta_hat_x, bxse = sqrt(D$sigma2_x), by = D$beta_hat_y, byse = sqrt(D$sigma2_y))
D$weights <- D$sigma2_y^-1
D$weights2 <- (D$beta_hat_y^2/D$beta_hat_x^2*D$sigma2_x + D$sigma2_y)^-1

# Plot least-squares and quantile regression line
plot(D$beta_hat_x, 
     D$beta_hat_y, 
     col = sapply(1:length(D$beta_hat_x), function(j) ifelse(j %in% D$valid.IVs, "black", "red")), 
     pch = 19
    )
abline(a = 0, b = 0.2, col = "grey", lwd = 2, lty = 2)
abline(lm(D$beta_hat_y ~ -1 + D$beta_hat_x, weights = D$weights), col = "orange", lwd = 2, lty = 2)
abline(lm(D$beta_hat_y_Egger ~ D$beta_hat_x_Egger, weights = D$weights), col = "darkgreen", lwd = 2, lty = 2)
abline(rq(D$beta_hat_y ~ -1 + D$beta_hat_x, weights = D$weights, tau = 0.3), col = "purple", lwd = 2, lty = 2)
abline(rq(D$beta_hat_y ~ -1 + D$beta_hat_x, weights = D$weights, tau = 0.5), col = "blue", lwd = 2, lty = 2)
abline(rq(D$beta_hat_y ~ -1 + D$beta_hat_x, weights = D$weights, tau = 0.7), col = "brown", lwd = 2, lty = 2)

# Add legend
legend("topleft", legend = c("y = theta * x", "IVW", "Egger", "QR tau=0.3", "QR tau=0.5", "QR tau=0.7"),
       col = c("grey", "orange", "darkgreen", "purple", "blue", "brown"), lwd = 2, lty = 2)

# # Calculate Cochran's Q statistic
# Q <- function(D, IVs = NULL, tau){
#   if (is.null(IVs)) IVs <- 1:length(D$beta_hat_y)
#   if (length(IVs) < 2) stop("Not enough valid IVs")
#   
#   # Calculate Cochran's Q statistic
#   df = length(IVs) - 1
#   qr.fit <- summary(rq(D$beta_hat_y[IVs] ~ -1 + D$beta_hat_x[IVs], weights = D$weights[IVs], tau = tau))
#   Q <- (D$beta_hat_y[IVs]/D$beta_hat_x[IVs] - qr.fit$coefficients)^2 %*% (D$weights[IVs] * D$beta_hat_x[IVs]^2)
#   RSE <- sqrt(Q/df)
#   return(list(RSE = c(RSE), df = df))
# }

# Fused lasso
fused_lasso.fit <- mr_fused_lasso(D$beta_hat_y, D$beta_hat_x, D$sigma2_y, nlambda = 10) 

# Extract fused estimates along the penalization path
thetas.df <- data.frame(t(lasso.fit$beta)) %>% round(4)
colnames(thetas.df) <- paste0("IV", 1:n)

# Compute RSS
RSS <- colSums((lasso.fit$y - lasso.fit$fit)^2)
AIC <- n * log(RSS / n) + 2 * lasso.fit$df
BIC <- n * log(RSS / n) + lasso.fit$df * log(n)

# Observed ratios
ratios_df <- data.frame(theta = sort(r)) %>%
  mutate(IV = paste0("IV", 1:n)) %>%
  mutate(IV = factor(IV, levels = paste0("IV", 1:n)))

# Plot distribution of estimated thetas with real values (in blue)
thetas.df %>% 
  pivot_longer(cols = everything(), names_to = "IV", values_to = "theta") %>%
  mutate(IV = factor(IV, levels = paste0("IV", 1:n))) %>%
  ggplot(aes(x = IV, y=theta)) +
  geom_boxplot() +
  geom_point(data = ratios_df, aes(x = IV, y = theta), 
             color = "blue", size = 2.5, shape = 18) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(y = "theta")

# Plot estimated thetas along the lasso path
path.df <- thetas.df %>% 
    pivot_longer(cols = everything(), names_to = "IV", values_to = "theta") %>%
    mutate(IV = factor(IV, levels = paste0("IV", 1:n))) %>%
    mutate(lambda = rep(paste0("lambda", 1:nrow(thetas.df)), each = n)) %>%
    mutate(lambda = factor(lambda, levels = (paste0("lambda", nrow(thetas.df):1)))) 

path.df %>%
  count(lambda, theta) %>%
  ggplot(aes(x = lambda, y=theta)) +
  geom_line(data=path.df, aes(group = IV), color = "grey", lty = 1) +
  geom_point(aes(size = n)) +
  scale_size_continuous(
    name = "Cluster size",    # Legend title
    range = c(1, 6),         # Size of points (min to max)
    breaks = c(1, 5, 10, 20, 25, 30), # Values shown in legend
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(y = expression(hat(theta)[lambda]), x = "Generalized lasso path")

# Plot least-squares and quantile regression line
theta_hat <- unlist(thetas.df[which.min(AIC), ])
D_ <- list()
D_$valid.IVs <- match(D$valid.IVs, order(r))
D_$beta_hat_x <- D$beta_hat_x[order(r)]
D_$beta_hat_y <- D_$beta_hat_x * theta_hat
D_$weights <- D$weights[order(r)]

plot(D_$beta_hat_x , 
     D_$beta_hat_y, 
     col = sapply(1:length(D$beta_hat_x), function(j) ifelse(j %in% D$valid.IVs, "black", "red"))[order(r)], 
     pch = 19
)
abline(a = 0, b = 0.2, col = "blue", lwd = 2, lty = 2)
for (i in 1:length(unique(theta_hat))){
  abline(a = 0, b = unique(theta_hat)[i], col = "grey", lwd = 2, lty = 2)
}

# Calculate median estimate using the biggest cluster
theta_mode <- names(sort(table(theta_hat), decreasing = TRUE)[1]) %>% as.numeric()
IVs <- order(r)[which(theta_hat == theta_mode)]
abline(rq(D$beta_hat_y[IVs] ~ -1 + D$beta_hat_x[IVs], weights = D$weights[IVs], tau = 0.5), col = "orange", lwd = 2, lty = 2)
