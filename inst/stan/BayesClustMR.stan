data {
  int<lower=1> N;                    // number of observations
  int<lower=1> K;                    // number of mixture components
  vector[N] theta;                   // observed ratios
  vector<lower=0>[N] b;              // known scale parameters for each observation
}

parameters {
  simplex[K] pi;                  // mixing proportions
  ordered[K] mu;                  // component locations (ordered for identifiability)
}

model {
  vector[K] log_pi = log(pi);  // cache log calculation
  
  // Priors
  mu ~ normal(0, 2);
  
  // Likelihood with observation-specific scale parameters
  for (n in 1:N) {
    vector[K] lps = log_pi;
    for (k in 1:K) {
      lps[k] += double_exponential_lpdf(theta[n] | mu[k], b[n]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  vector[N] log_lik;     // log-likelihood for each observation
  matrix[N, K] r;        // responsibilities
  
  // Calculate log-likelihood and assign components
  for (n in 1:N) {
    vector[K] lps = log(pi);
    for (k in 1:K) {
      lps[k] += double_exponential_lpdf(theta[n] | mu[k], b[n]);
    }
    log_lik[n] = log_sum_exp(lps);
    r[n, ] = to_row_vector(softmax(lps));
  }
}
