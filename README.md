BayesFuseMR
================

### Overview

This repository provides the **BayesFuseMR** R package to implement our
Bayesian fused shrinkage model for MR studies, along with scripts to
reproduce the simulations and results described in our ongoing
manuscript.

### Installation

You first need to install the latest release of the **cmdstanr** R
package with

``` r
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

If you don’t already have CmdStan installed then, in addition to
installing the R package, it is also necessary to install CmdStan using
CmdStanR’s `install_cmdstan()` function. A suitable C++ toolchain is
also required. Instructions are provided in the [Getting started with
CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) vignette.
You can install the **BayesFuseMR** R package from GitHub using the
`devtools` package:

``` r
devtools::install_github("julstpierre/BayesFuseMR")
```

We can now load the package like any other R package. We’ll also load
the **bayesplot** and **dplyr** packages to use later in examples.

``` r
library(BayesFuseMR)
library(bayesplot)
library(dplyr)
```

### Load example data

We provide a dummy dataset in the package to illustrate how to use the
`fit_BayesFuseMR()` function. This dataset contains simulated summary
statistics from 30 SNPs for one continuous exposure and one continuous
outcome.

``` r
data("dummy_data")
head(dummy_data)
```

    ##     beta_hat_y beta_hat_x    sigma_y    sigma_x
    ## 1  0.030079854 0.14838236 0.02355933 0.02218627
    ## 2  0.001951586 0.18661405 0.02674638 0.02513565
    ## 3  0.285058236 0.20056546 0.02909100 0.02766173
    ## 4  0.046815342 0.14481593 0.02278348 0.02141118
    ## 5  0.043159749 0.20483177 0.03223177 0.03071550
    ## 6 -0.021999053 0.04590455 0.02459339 0.02313552

### Running MCMC

The function `fit_BayesFuseMR()` fits a Bayesian fused shrinkage model
using MCMC. It takes as input the estimated effect sizes for the
exposure and outcome with their standard errors, and returns a fitted
`cmdstanr` object. We can use a seed when sampling for reproducibility.

``` r
fit <- fit_BayesFuseMR(
  beta_hat_y = dummy_data$beta_hat_y,
  beta_hat_x = dummy_data$beta_hat_x,
  sigma_y = dummy_data$sigma_y,
  sigma_x = dummy_data$sigma_x,
  snps = paste0("SNP", 1:nrow(dummy_data)),
  seed = 123
)
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 4.7 seconds.
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 4.7 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 4.9 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 5.3 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 4.9 seconds.
    ## Total execution time: 5.4 seconds.

### Posterior summary statistics

You can use the `summary` method to get posterior summaries of all
parameters. Here, we combine the SNP identifiers with the posterior
summaries of the `theta` parameter, which represents the ratio estimate
for each variant.

``` r
fit_summary <- cbind(snp = fit$snps, fit$posterior$summary(variables = c("theta")))
head(fit_summary)
```

    ##     snp variable        mean     median        sd       mad         q5
    ## 1 SNP23 theta[1] -0.44496524 -0.4425505 0.1140987 0.1129497 -0.6377542
    ## 2 SNP10 theta[2] -0.33773368 -0.3394260 0.1104323 0.1102573 -0.5247162
    ## 3 SNP26 theta[3] -0.28713810 -0.2882315 0.1439096 0.1389515 -0.5199317
    ## 4  SNP6 theta[4] -0.22774694 -0.2277370 0.1576992 0.1490784 -0.4935032
    ## 5 SNP29 theta[5] -0.14473035 -0.1443235 0.1487275 0.1468406 -0.3837784
    ## 6 SNP13 theta[6] -0.06618775 -0.0654981 0.1186513 0.1198953 -0.2632873
    ##           q95      rhat ess_bulk ess_tail
    ## 1 -0.26168190 1.0000732 4311.110 3709.090
    ## 2 -0.15724660 1.0003949 3484.147 3729.837
    ## 3 -0.04666115 0.9998433 3673.113 3041.839
    ## 4  0.03035697 1.0009753 3241.509 3262.469
    ## 5  0.09582668 1.0010242 3093.753 3034.730
    ## 6  0.12762665 1.0000261 5100.352 3976.570

### Posterior draws

You can extract the posterior draws for further analysis or
visualization. The `draws()` method returns a `draws_df` object, which
can be converted to a data frame for easier manipulation. For example,
we can look at the trace plots and density plots of the posterior draws
for the `theta` parameter to ensure convergence of the MCMC chains.

``` r
# Trace plots for theta
mcmc_trace(fit$posterior$draws("theta"))
```

![](README_files/figure-gfm/trace%20plots-1.png)<!-- -->

``` r
# Density plots for theta
mcmc_dens_chains(fit$posterior$draws("theta"))
```

![](README_files/figure-gfm/density%20plots-1.png)<!-- -->

### Clustering

The function `clust_BayesFuseMR()` identifies clusters and computes
weighted cluster means, standard errors, and confidence intervals. We
can select the number of clusters using the FDR argument which controls
the false discovery rate threshold for cluster selection. Here, we set
it to 0.4.

``` r
clust_BayesFuseMR(fit, FDR = 0.4) %>% 
  print(n=Inf)
```

    ## # A tibble: 30 × 10
    ##    snp   cluster cluster_size cluster_mean cluster_se cluster_lcl cluster_ucl
    ##    <chr>   <int>        <int>        <dbl>      <dbl>       <dbl>       <dbl>
    ##  1 SNP23       1            1       -0.445     0.114      -0.669      -0.221 
    ##  2 SNP10       2            2       -0.319     0.0876     -0.491      -0.147 
    ##  3 SNP26       2            2       -0.319     0.0876     -0.491      -0.147 
    ##  4 SNP6        3            1       -0.228     0.158      -0.537       0.0813
    ##  5 SNP29       4            1       -0.145     0.149      -0.436       0.147 
    ##  6 SNP13       5           19        0.151     0.0262      0.0992      0.202 
    ##  7 SNP9        5           19        0.151     0.0262      0.0992      0.202 
    ##  8 SNP2        5           19        0.151     0.0262      0.0992      0.202 
    ##  9 SNP12       5           19        0.151     0.0262      0.0992      0.202 
    ## 10 SNP8        5           19        0.151     0.0262      0.0992      0.202 
    ## 11 SNP28       5           19        0.151     0.0262      0.0992      0.202 
    ## 12 SNP16       5           19        0.151     0.0262      0.0992      0.202 
    ## 13 SNP7        5           19        0.151     0.0262      0.0992      0.202 
    ## 14 SNP18       5           19        0.151     0.0262      0.0992      0.202 
    ## 15 SNP20       5           19        0.151     0.0262      0.0992      0.202 
    ## 16 SNP15       5           19        0.151     0.0262      0.0992      0.202 
    ## 17 SNP14       5           19        0.151     0.0262      0.0992      0.202 
    ## 18 SNP25       5           19        0.151     0.0262      0.0992      0.202 
    ## 19 SNP1        5           19        0.151     0.0262      0.0992      0.202 
    ## 20 SNP5        5           19        0.151     0.0262      0.0992      0.202 
    ## 21 SNP22       5           19        0.151     0.0262      0.0992      0.202 
    ## 22 SNP11       5           19        0.151     0.0262      0.0992      0.202 
    ## 23 SNP4        5           19        0.151     0.0262      0.0992      0.202 
    ## 24 SNP21       5           19        0.151     0.0262      0.0992      0.202 
    ## 25 SNP19       6            2        0.465     0.0792      0.310       0.621 
    ## 26 SNP24       6            2        0.465     0.0792      0.310       0.621 
    ## 27 SNP17       7            1        0.526     0.0943      0.341       0.711 
    ## 28 SNP27       8            1        0.699     0.101       0.502       0.896 
    ## 29 SNP30       9            1        1.10      0.197       0.716       1.49  
    ## 30 SNP3       10            1        1.25      0.155       0.950       1.56  
    ## # ℹ 3 more variables: theta_mean <dbl>, theta_se <dbl>, theta_pep <dbl>
