library(ggplot2)
library(dplyr)

# Prepare data for the lines and custom legend labels
dat <- cbind(beta_hat_y = D$beta_hat_y, 
             beta_hat_y_se = sqrt(D$sigma2_y), 
             beta_hat_x = D$beta_hat_x, 
             beta_hat_x_se = sqrt(D$sigma2_x), 
             res_em$results
             )

cluster_lines <- dat %>%
  distinct(cluster, cluster_quantile) %>%
  mutate(
    cluster_label = paste0("Classe ", cluster, " (", 
                           round(cluster_quantile, 3), ")")
  )

# Merge labels back to main data for consistent coloring
dat <- dat %>%
  left_join(cluster_lines %>% select(cluster, cluster_label), by = "cluster")

ggplot(dat, aes(x = beta_hat_x, y = beta_hat_y)) +
  # Add error bars (y)
  geom_errorbar(aes(ymin = beta_hat_y - beta_hat_y_se,
                    ymax = beta_hat_y + beta_hat_y_se,
                    color = factor(cluster_label)),
                width = 0.002, alpha = 0.6) +
  # Add error bars (x)
  geom_errorbarh(aes(xmin = beta_hat_x - beta_hat_x_se,
                     xmax = beta_hat_x + beta_hat_x_se,
                     color = factor(cluster_label)),
                 height = 0.002, alpha = 0.6) +
  # Add points with size proportional to inclusion probability
  geom_point(aes(size = probability, color = factor(cluster_label)), shape = 1,
             alpha = 0.8) +
  # Add cluster lines with slope = cluster_quantile
  geom_abline(aes(slope = cluster_quantile, intercept = 0,
                  color = factor(cluster_label)),
              linetype = "dashed", linewidth = 0.5) +
  scale_size_continuous(name = "Probabilité d'inclusion") +
  scale_color_discrete(name = "Classe (effet causal médian)") +
  labs(
    x = "Effets génétiques sur la variable d'exposition",
    y = "Effets génétiques sur la variable réponse",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12)
  )
