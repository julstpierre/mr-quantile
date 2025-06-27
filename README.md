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
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
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
fit <- fit_BayesFuseMR(beta_hat_y = dummy_data$beta_hat_y, beta_hat_x = dummy_data$beta_hat_x,
    sigma_y = dummy_data$sigma_y, sigma_x = dummy_data$sigma_x, snps = paste0("SNP",
        1:nrow(dummy_data)))
```

    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 4.4 seconds.
    ## Chain 4 finished in 4.4 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 4.5 seconds.
    ## Chain 3 finished in 4.4 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 4.4 seconds.
    ## Total execution time: 4.6 seconds.

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
    ## 1 SNP23 theta[1] -0.44200066 -0.43902500 0.1135344 0.1113040 -0.6331745
    ## 2 SNP10 theta[2] -0.33722183 -0.33712200 0.1127035 0.1106516 -0.5217042
    ## 3 SNP26 theta[3] -0.28910864 -0.29193250 0.1463468 0.1417655 -0.5224859
    ## 4  SNP6 theta[4] -0.22556145 -0.22294300 0.1583803 0.1516396 -0.4831539
    ## 5 SNP29 theta[5] -0.14607205 -0.14528000 0.1485054 0.1435918 -0.3917273
    ## 6 SNP13 theta[6] -0.06556177 -0.06421955 0.1193351 0.1166066 -0.2632553
    ##           q95     rhat ess_bulk ess_tail
    ## 1 -0.25971925 1.000226 3814.300 3798.541
    ## 2 -0.15131970 1.000357 4055.049 3760.828
    ## 3 -0.04457441 1.000593 3913.688 3406.526
    ## 4  0.02406007 1.000998 3417.889 3232.759
    ## 5  0.08739021 1.001012 3875.691 3325.002
    ## 6  0.12733650 1.000809 4718.847 3781.367

The function `clust_BayesFuseMR()` identifies clusters and computes
weighted cluster means, standard errors, and confidence intervals.

``` r
clust_BayesFuseMR(fit, FDR = 0.3)
```

    ## # A tibble: 30 × 10
    ##    snp   cluster cluster_size cluster_mean cluster_se cluster_lcl cluster_ucl
    ##    <chr>   <int>        <int>        <dbl>      <dbl>       <dbl>       <dbl>
    ##  1 SNP23       1            1       -0.442     0.114      -0.665       -0.219
    ##  2 SNP10       2            3       -0.297     0.0778     -0.449       -0.144
    ##  3 SNP26       2            3       -0.297     0.0778     -0.449       -0.144
    ##  4 SNP6        2            3       -0.297     0.0778     -0.449       -0.144
    ##  5 SNP29       3           20        0.145     0.0257      0.0942       0.195
    ##  6 SNP13       3           20        0.145     0.0257      0.0942       0.195
    ##  7 SNP9        3           20        0.145     0.0257      0.0942       0.195
    ##  8 SNP2        3           20        0.145     0.0257      0.0942       0.195
    ##  9 SNP12       3           20        0.145     0.0257      0.0942       0.195
    ## 10 SNP8        3           20        0.145     0.0257      0.0942       0.195
    ## # ℹ 20 more rows
    ## # ℹ 3 more variables: theta_mean <dbl>, theta_se <dbl>, theta_pep <dbl>
