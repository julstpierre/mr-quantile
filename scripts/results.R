# Load libraries
library(dplyr)
library(data.table)
library(tidyr)

# Load results
res_theta <- rbind(fread("results/theta_scenario1.txt", sep = "\t"),
               fread("results/theta_scenario2.txt", sep = "\t"),
               fread("results/theta_scenario3.txt", sep = "\t")
             ) %>% 
             mutate(method = gsub("theta_", "", V1)) %>%
             select(-V1)

# Table 1 - theta = 0.2
filter(res_theta, theta == 0.2) %>%
  group_by(scenario, frac, method) %>%
summarise(Mean = mean(Estimate, na.rm = TRUE), SD = sd(Estimate, na.rm = TRUE), SE = mean(`Std. Error`, na.rm = TRUE), Power = mean(`Pr(>|t|)` < 0.05, na.rm = TRUE)) %>%
  pivot_wider(names_from = frac, values_from = c(Mean, SD, SE, Power), names_vary = "slowest") %>%
  arrange(scenario, match(method, c('IVW', 'MR_IVW', 'Egger', 'MR_Egger','Median', 'lasso')))

# Table 2 - theta = 0
filter(res_theta, theta == 0) %>%
  group_by(scenario, frac, method) %>%
  summarise(Mean = mean(Estimate, na.rm = TRUE), SD = sd(Estimate, na.rm = TRUE), SE = mean(`Std. Error`, na.rm = TRUE), `Type I` = mean(`Pr(>|t|)` < 0.05, na.rm = TRUE)) %>%
  pivot_wider(names_from = frac, values_from = c(Mean, SD, SE, `Type I`), names_vary = "slowest") %>%
  arrange(scenario, match(method, c('IVW', 'MR_IVW','Egger', 'MR_Egger','Median', 'lasso')))
  