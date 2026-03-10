source("utils.R")

library(dplyr)
library(purrr)
library(TwoSampleMR)
library(MendelianRandomization)
library(MRcML)
library(tidyr)

methods <- list(
  mr_ivw       = fit_mr_ivw,
  mr_median    = fit_mr_median,
  mr_egger     = fit_mr_egger,
  mr_mode      = fit_mr_mode,
  mr_raps      = fit_mr_raps,
  mr_lasso     = fit_mr_lasso,
  mr_contmix   = fit_mr_contmix,
  mr_mix       = fit_mr_mix,
  mr_quantreg  = fit_mr_quantreg,
  mr_cML       = fit_mr_cML,
  mr_presso    = fit_mr_presso,
  mr_cML_DP    = fit_mr_cML_DP
)

method_labels <- c(
  mr_cML        = "cML-MA-BIC",
  mr_cML_DP     = "cML-MA-BIC-DP",
  mr_contmix    = "MR-ContMix",
  mr_egger      = "MR-Egger",
  mr_ivw        = "MR-IVW",
  mr_lasso      = "MR-Lasso",
  mr_median     = "MR-Weighted-Median",
  mr_mix        = "MR-Mix",
  mr_mode       = "MR-Weighted-Mode",
  mr_presso     = "MR-PRESSO",
  mr_quantreg   = "MR-Quantile",
  mr_raps       = "MR-RAPS"
)

D <- D_S41[1:100]
res30 <- imap_dfr(methods, function(f, method) {
  map_dbl(D, ~ system.time(f(.x))["elapsed"]) |>
    tibble(
      method = method,
      time   = _
    )
}) |>
  mutate(p = 30)

D <- D_S42[1:100]
res100 <- imap_dfr(methods, function(f, method) {
  map_dbl(D, ~ system.time(f(.x))["elapsed"]) |>
    tibble(
      method = method,
      time   = _
    )
}) |>
  mutate(p = 100)

# Print results
fmt_mean_sd <- function(mean, sd, thresh = 0.01) {
  if (mean < thresh && sd < thresh) {
    sprintf("%.3f (%.3f)", mean, sd)
  } else if (mean < thresh) {
    sprintf("%.3f (%.2f)", mean, sd)
  } else  if (sd < thresh) {
    sprintf("%.2f (%.3f)", mean, sd)
  } else {
    sprintf("%.2f (%.2f)", mean, sd)
  }
}

bind_rows(res30, res100) |>
  group_by(p, method) |>
  summarise(
    mean   = mean(time),
    sd     = sd(time),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = mapply(fmt_mean_sd, mean, sd)) %>%
  transmute(method, p, mean_sd) %>%
  mutate(
    method = factor(method,
                    levels = names(method_labels),
                    labels = method_labels)
  ) %>%
  pivot_wider(names_from = p, values_from = mean_sd) %>%
  arrange(`100`)
