// Provides the STAN implementation of the MR-Horse method described in 
// Grant, AJ and Burgess, S (2024). A Bayesian approach to Mendelian randomization 
// using summary statistics in the univariable and multivariable settings with correlated pleiotropy. 
// American Journal of Human Genetics. 111(1):165-180. doi: https://doi.org/10.1016/j.ajhg.2023.12.002

data {
  int<lower=0> p;              // number of instruments
  vector[p] Gamma_hat_y;       // SNP-outcome associations
  vector[p] gamma_hat_x;       // SNP-exposure associations
  vector<lower=0>[p] sigma_y;  // standard errors of Gamma_hat_y
  vector<lower=0>[p] sigma_x;  // standard errors of gamma_hat_x
  real<lower=0> scale_global;  // scale for the half cauchy prior for tau
  real<lower=0> slab_scale;    // slab scale for the regularized horseshoe
  real<lower=0> slab_df;       // slab degrees of freedom for the regularized horseshoe
}

parameters {
  real beta;                     // causal effect                
  vector[2*p] z;                 // auxiliary variable
  real mu;                       // mean of gamma_x          
  real<lower=0> sx;              // standard deviation of gamma_x
  vector<lower=0, upper=1>[p] r; // correlation hyperparameter
  vector<lower=0>[p] aux1_local;
  vector<lower=0>[p] aux2_local;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  real<lower=0> caux ;
}

transformed parameters{
  vector[p] gamma_x; // unpenalized regression parameters
  vector[p] alpha; // penalized regression parameters
  vector<lower=-1, upper=1>[p] rho; // correlation parameters
  vector<lower=0>[p] lambda; // local scale parameters
  vector<lower=0>[p] lambda_tilde; // ’ truncated ’ local shrinkage parameter
  real<lower=0> tau; // global scale parameter
  real<lower=0> c; // slab scale
  vector[p] theta; // linear predictor
  
  // The order of the statements is important, for example defining alpha before gamma_x ! 
  rho = 2 * r - 1;
  lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global;
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt( c^2 * square(lambda) ./ ( c^2 + tau^2 * square(lambda) ) );
  
  alpha = z[p+1:] .* lambda_tilde * tau ;
  gamma_x = mu + rho * sx ./ ( tau * lambda_tilde ) .* alpha + z[1:p] .* sqrt( 1.0 - square(rho) ) * sx;
  theta = beta * gamma_x + alpha;
}

model {
  // Normal likelihood
  for (i in 1:p){
    Gamma_hat_y[i] ~ normal(theta[i] * gamma_hat_x[i], sigma_y[i]);
  }
  
  // Priors for hyperparameters
  z ~ normal(0, 1);
  mu ~ normal(0, 1);
  sx ~ normal(0, 1);
  beta ~ normal(0, 1);
  r ~ beta(10, 10);
  
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5, 0.5);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5, 0.5);
  caux ~ inv_gamma(0.5* slab_df , 0.5* slab_df);
}
