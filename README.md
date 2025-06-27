### Overview

This repository provides the **BayesFuseMR** R package to implement our
Bayesian fused shrinkage model for MR studies, along with scripts to
reproduce the simulations and results described in our ongoing
manuscript.

### Installation

You first need to install the latest release of the **cmdstanr** R
package with

    install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

If you don’t already have CmdStan installed then, in addition to
installing the R package, it is also necessary to install CmdStan using
CmdStanR’s `install_cmdstan()` function. A suitable C++ toolchain is
also required. Instructions are provided in the [Getting started with
CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) vignette.
You can install the **BayesFuseMR** R package from GitHub using the
`devtools` package:

    devtools::install_github("julstpierre/BayesFuseMR")

We can now load the package like any other R package. We’ll also load
the **bayesplot** and **posterior** packages to use later in examples.

    library(BayesFuseMR)
    library(bayesplot)
    library(posterior)

### Load example data

We provide a dummy dataset in the package to illustrate how to use the
`fit_BayesFuseMR()` function. This dataset contains simulated summary
statistics from p=30 SNPs for one continuous exposure and one continuous
outcome.

    data("dummy_data")
    head(dummy_data)
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

    fit <- fit_BayesFuseMR(beta_hat_y = dummy_data$beta_hat_y, beta_hat_x = dummy_data$beta_hat_x,
        sigma_y = dummy_data$sigma_y, sigma_x = dummy_data$sigma_x, snps = paste0("SNP", 1:nrow(dummy_data)))
    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 4.5 seconds.
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 4.7 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    ## Chain 1 finished in 4.8 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 6.0 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 5.0 seconds.
    ## Total execution time: 6.1 seconds.

### Posterior summaries

You can use the `summary` method to get posterior summaries of all
parameters. Here, we combine the SNP identifiers with the posterior
summaries of the `theta` parameter, which represents the ratio estimate
for each variant.

    fit_summary <- cbind(snp = fit$snps, fit$posterior$summary(variables = c("theta")))
    head(fit_summary)
    ##     snp variable        mean      median        sd       mad         q5         q95
    ## 1 SNP23 theta[1] -0.44118842 -0.43902150 0.1145044 0.1172195 -0.6298249 -0.25578430
    ## 2 SNP10 theta[2] -0.33777858 -0.33716050 0.1124481 0.1168897 -0.5193524 -0.15221890
    ## 3 SNP26 theta[3] -0.28879679 -0.28800700 0.1434637 0.1365741 -0.5289193 -0.05810620
    ## 4  SNP6 theta[4] -0.22706267 -0.22349750 0.1627209 0.1562001 -0.4960164  0.04189190
    ## 5 SNP29 theta[5] -0.14706147 -0.15136200 0.1511152 0.1469348 -0.3904361  0.09751309
    ## 6 SNP13 theta[6] -0.06764401 -0.06952295 0.1206429 0.1187252 -0.2646341  0.13269300
    ##        rhat ess_bulk ess_tail
    ## 1 1.0012190 4134.726 3830.906
    ## 2 0.9995740 4241.325 3786.749
    ## 3 0.9998649 4740.102 3616.850
    ## 4 1.0000074 4955.907 3634.807
    ## 5 1.0011181 4635.561 3760.392
    ## 6 1.0000804 5775.384 3447.917

The function `clust_BayesFuseMR()` identifies clusters and computes
weighted cluster means, standard errors, and confidence intervals.

    clust_BayesFuseMR(fit, FDR = 0.3)
    ## # A tibble: 30 × 10
    ##    snp   cluster cluster_size cluster_mean cluster_se cluster_lcl cluster_ucl theta_mean
    ##    <chr>   <int>        <int>        <dbl>      <dbl>       <dbl>       <dbl>      <dbl>
    ##  1 SNP23       1            1       -0.441     0.115      -0.666       -0.217   -0.441  
    ##  2 SNP10       2            4       -0.267     0.0691     -0.402       -0.131   -0.338  
    ##  3 SNP26       2            4       -0.267     0.0691     -0.402       -0.131   -0.289  
    ##  4 SNP6        2            4       -0.267     0.0691     -0.402       -0.131   -0.227  
    ##  5 SNP29       2            4       -0.267     0.0691     -0.402       -0.131   -0.147  
    ##  6 SNP13       3           19        0.150     0.0266      0.0978       0.202   -0.0676 
    ##  7 SNP9        3           19        0.150     0.0266      0.0978       0.202   -0.0271 
    ##  8 SNP2        3           19        0.150     0.0266      0.0978       0.202    0.00335
    ##  9 SNP12       3           19        0.150     0.0266      0.0978       0.202    0.0342 
    ## 10 SNP8        3           19        0.150     0.0266      0.0978       0.202    0.0601 
    ## # ℹ 20 more rows
    ## # ℹ 2 more variables: theta_se <dbl>, theta_pep <dbl>
