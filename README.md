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
##     beta_hat_y beta_hat_x    sigma_y    sigma_x
## 1  0.030079854 0.14838236 0.02355933 0.02218627
## 2  0.001951586 0.18661405 0.02674638 0.02513565
## 3  0.285058236 0.20056546 0.02909100 0.02766173
## 4  0.046815342 0.14481593 0.02278348 0.02141118
## 5  0.043159749 0.20483177 0.03223177 0.03071550
## 6 -0.021999053 0.04590455 0.02459339 0.02313552
```

### Fit the model

The function `fit_BayesFuseMR()` fits a Bayesian fused shrinkage model
using MCMC. It takes as input the estimated effect sizes for the
exposure and outcome with their standard errors, and returns a fitted
`cmdstanr` object.

``` r
fit <- fit_BayesFuseMR(beta_hat_y = dummy_data$beta_hat_y, beta_hat_x = dummy_data$beta_hat_x,
    sigma_y = dummy_data$sigma_y, sigma_x = dummy_data$sigma_x, snps = paste0("SNP",
        1:nrow(dummy_data)))
## Running MCMC with 4 parallel chains...
## 
## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
## Chain 2 finished in 4.4 seconds.
## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
## Chain 3 finished in 4.6 seconds.
## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
## Chain 1 finished in 5.3 seconds.
## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
## Chain 4 finished in 6.1 seconds.
## 
## All 4 chains finished successfully.
## Mean chain execution time: 5.1 seconds.
## Total execution time: 6.2 seconds.
```

### Posterior summaries

You can use the `summary` method to get posterior summaries of all
parameters. Here, we combine the SNP identifiers with the posterior
summaries of the `theta` parameter, which represents the ratio estimate
for each variant.

``` r
fit_summary <- cbind(snp = fit$snps, fit$posterior$summary(variables = c("theta")))
head(fit_summary)
##     snp variable        mean     median        sd       mad         q5
## 1 SNP23 theta[1] -0.44127635 -0.4395105 0.1140586 0.1112995 -0.6352581
## 2 SNP10 theta[2] -0.33818467 -0.3387580 0.1106949 0.1076909 -0.5153921
## 3 SNP26 theta[3] -0.28617316 -0.2869860 0.1449396 0.1395905 -0.5258101
## 4  SNP6 theta[4] -0.22731696 -0.2281710 0.1598894 0.1538724 -0.4817818
## 5 SNP29 theta[5] -0.14720842 -0.1468590 0.1534328 0.1475417 -0.3995317
## 6 SNP13 theta[6] -0.06783503 -0.0663992 0.1228092 0.1199670 -0.2722240
##           q95      rhat ess_bulk ess_tail
## 1 -0.25749360 0.9999044 3586.668 3674.456
## 2 -0.15152565 1.0003056 3312.986 3550.556
## 3 -0.04499278 1.0011981 3634.540 3103.999
## 4  0.02961725 1.0016566 3397.597 3345.953
## 5  0.10018100 1.0000568 3574.701 3598.555
## 6  0.13620760 1.0008568 5006.497 3905.777
```

The function `clust_BayesFuseMR()` identifies clusters and computes
weighted cluster means, standard errors, and confidence intervals.

``` r
clust_BayesFuseMR(fit, FDR = 0.3)
## # A tibble: 30 × 10
##    snp   cluster cluster_size cluster_mean cluster_se cluster_lcl cluster_ucl
##    <chr>   <int>        <int>        <dbl>      <dbl>       <dbl>       <dbl>
##  1 SNP23       1            1       -0.441     0.114      -0.665       -0.218
##  2 SNP10       2            3       -0.298     0.0771     -0.449       -0.147
##  3 SNP26       2            3       -0.298     0.0771     -0.449       -0.147
##  4 SNP6        2            3       -0.298     0.0771     -0.449       -0.147
##  5 SNP29       3           20        0.143     0.0260      0.0925       0.194
##  6 SNP13       3           20        0.143     0.0260      0.0925       0.194
##  7 SNP9        3           20        0.143     0.0260      0.0925       0.194
##  8 SNP2        3           20        0.143     0.0260      0.0925       0.194
##  9 SNP12       3           20        0.143     0.0260      0.0925       0.194
## 10 SNP8        3           20        0.143     0.0260      0.0925       0.194
## # ℹ 20 more rows
## # ℹ 3 more variables: theta_mean <dbl>, theta_se <dbl>, theta_pep <dbl>
```
