MRQuantile
================

### Overview

This repository provides the **MRQuantile** R package to implement our
weighted quantile regression model for two samples MR studies, along
with scripts to reproduce the simulations and results described in our
ongoing manuscript.

### Installation

You can install the **MRQuantile** R package from GitHub using the
`devtools` package:

``` r
devtools::install_github("julstpierre/MRQuantile")
```

We can now load the package like any other R package. We’ll also load
the **dplyr** package to use later in examples.

``` r
library(MRQuantile)
library(dplyr)
```

### Load example data

We provide a dummy dataset in the package to illustrate how to use the
`mr_quantile()` function. This dataset contains simulated summary
statistics from 30 SNPs for one continuous exposure and one continuous
outcome.

``` r
head(dummy_data)
```

    ##     beta_hat_y beta_hat_x    sigma_y    sigma_x
    ## 1  0.030079854 0.14838236 0.02355933 0.02218627
    ## 2  0.001951586 0.18661405 0.02674638 0.02513565
    ## 3  0.285058236 0.20056546 0.02909100 0.02766173
    ## 4  0.046815342 0.14481593 0.02278348 0.02141118
    ## 5  0.043159749 0.20483177 0.03223177 0.03071550
    ## 6 -0.021999053 0.04590455 0.02459339 0.02313552

### Fitting the model

The function `mr_quantile()` fits a weighted quantile regression model
using the asymmetric Laplace distribution. It takes as input the
estimated effect sizes for the exposure and outcome with their
variances, and returns a list with `coefficients`.

``` r
fit <- mr_quantile(
  beta_hat_y = dummy_data$beta_hat_y,
  beta_hat_x = dummy_data$beta_hat_x,
  sigma2_y = dummy_data$sigma_y^2,
  sigma2_x = dummy_data$sigma_x^2,
  verbose = TRUE
)
```

    ## Converged after 5 iterations
    ## Estimated tau: 0.4342752 
    ## Estimated lambda: 1.840049

``` r
fit$coefficients
```

    ##    Estimate   StdError         LCL       UCL
    ## 1 0.1284859 0.09033615 -0.04856965 0.3055415
