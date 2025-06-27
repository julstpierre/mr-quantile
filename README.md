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
the **bayesplot** and **posterior** packages to use later in examples.

``` r
library(BayesFuseMR)
library(bayesplot)
library(posterior)
```

### Load example data

We provide a dummy dataset in the package to illustrate how to use the
`fit_BayesFuseMR()` function. This dataset contains simulated summary
statistics from p=30 SNPs for one continuous exposure and one continuous
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

### Fit the model

The function `fit_BayesFuseMR()` fits a Bayesian fused shrinkage model
using MCMC. It takes as input the estimated effect sizes for the
exposure and outcome with their standard errors, and returns a fitted
`cmdstanr` object.

``` r
fit <- fit_BayesFuseMR(
  beta_hat_y = dummy_data$beta_hat_y,
  beta_hat_x = dummy_data$beta_hat_x,
  sigma_y = dummy_data$sigma_y,
  sigma_x = dummy_data$sigma_x,
  snps = paste0("SNP", 1:nrow(dummy_data))
)
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 4.7 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 4.9 seconds.
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 5.2 seconds.
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 5.9 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 5.2 seconds.
    ## Total execution time: 6.1 seconds.

### Posterior summaries

You can use the `summary` method to get posterior summaries of all
parameters. Here, we combine the SNP identifiers with the posterior
summaries of the `theta` parameter, which represents the ratio estimate
for each variant.

``` r
fit_summary <- cbind(snp = fit$snps, fit$posterior$summary(variables = c("theta")))
head(fit_summary)
```

    ##     snp variable        mean      median        sd       mad         q5
    ## 1 SNP23 theta[1] -0.44374277 -0.44172900 0.1100599 0.1095604 -0.6276431
    ## 2 SNP10 theta[2] -0.33805301 -0.33913900 0.1116392 0.1104159 -0.5241826
    ## 3 SNP26 theta[3] -0.28492620 -0.28706300 0.1456072 0.1405127 -0.5186931
    ## 4  SNP6 theta[4] -0.22759425 -0.22934800 0.1591705 0.1542756 -0.4865954
    ## 5 SNP29 theta[5] -0.14778729 -0.14526200 0.1500656 0.1479998 -0.3997963
    ## 6 SNP13 theta[6] -0.06890449 -0.06804025 0.1187223 0.1189872 -0.2651305
    ##           q95      rhat ess_bulk ess_tail
    ## 1 -0.26263985 1.0001693 4113.284 3584.727
    ## 2 -0.15445800 1.0003834 4315.544 3618.415
    ## 3 -0.04695237 1.0004753 3973.852 3329.009
    ## 4  0.03983069 0.9998311 4184.879 3687.291
    ## 5  0.09428248 1.0007180 4201.765 3223.944
    ## 6  0.12388585 0.9997957 5477.623 3325.693

The function `clust_BayesFuseMR()` identifies clusters and computes
weighted cluster means, standard errors, and confidence intervals.

``` r
clust_BayesFuseMR(fit, FDR = 0.3)
```

    ## # A tibble: 30 × 10
    ##    snp   cluster cluster_size cluster_mean cluster_se cluster_lcl cluster_ucl
    ##    <chr>   <int>        <int>        <dbl>      <dbl>       <dbl>       <dbl>
    ##  1 SNP23       1            1       -0.444     0.110      -0.659       -0.228
    ##  2 SNP10       2            3       -0.297     0.0774     -0.449       -0.145
    ##  3 SNP26       2            3       -0.297     0.0774     -0.449       -0.145
    ##  4 SNP6        2            3       -0.297     0.0774     -0.449       -0.145
    ##  5 SNP29       3           20        0.142     0.0258      0.0915       0.193
    ##  6 SNP13       3           20        0.142     0.0258      0.0915       0.193
    ##  7 SNP9        3           20        0.142     0.0258      0.0915       0.193
    ##  8 SNP2        3           20        0.142     0.0258      0.0915       0.193
    ##  9 SNP12       3           20        0.142     0.0258      0.0915       0.193
    ## 10 SNP8        3           20        0.142     0.0258      0.0915       0.193
    ## # ℹ 20 more rows
    ## # ℹ 3 more variables: theta_mean <dbl>, theta_se <dbl>, theta_pep <dbl>
