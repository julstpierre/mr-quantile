#----------------------------
# Mr-clust
#----------------------------

data.df <- dummy_data
data.df$ratios <- data.df$beta_hat_y / data.df$beta_hat_x
data.df$ratios.se <- data.df$sigma_y / abs(data.df$beta_hat_x)
ord <- order(data.df$ratios / data.df$ratios.se)
data.df <- data.df[ord, ]
res_em <- mr_clust_em(theta = data.df$ratios, 
                      theta_se = data.df$ratios.se, 
                      bx = data.df$beta_hat_x, 
                      by = data.df$beta_hat_y, 
                      bxse = data.df$sigma_x, 
                      byse = data.df$sigma_y,
                      obs_names = paste0("SNP", 1:p)
)
res_em$results$best
