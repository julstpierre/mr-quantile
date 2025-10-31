
data {
  int<lower=0> p;               // number of instruments
  vector[p] pi_scaled;          // Wald ratios (scaled)
  vector<lower=0>[p] sigma_pi;  // standard errors of ratios
  real<lower=0> scale_local;    // scale for the normal prior for theta[1]
  real<lower=0> scale_global;   // scale for the half cauchy prior for tau
  real<lower=0> slab_scale;     // slab scale for the regularized horseshoe
  real<lower=0> slab_df;        // slab degrees of freedom for the regularized horseshoe
}

parameters {
  real mu;                      // prior mean for theta[1]
  vector[p] z;                  // auxiliary variable
  vector<lower=0>[p-1] aux1_local;
  vector<lower=0>[p-1] aux2_local;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  real<lower=0> caux ;
}

transformed parameters{
  vector<lower=0>[p-1] lambda; // local scale parameters
  vector<lower=0>[p-1] lambda_tilde; // ’ truncated ’ local shrinkage parameter
  real<lower=0> tau; // global scale parameter
  real<lower=0> c; // slab scale
  vector[p] theta; // mean parameters
  vector[p-1] eta; // succesive differences
  
  // The order of the statements is important 
  lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global;
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt( c^2 * square(lambda) ./ ( c^2 + tau^2 * square(lambda) ) );
  
  theta[1] = mu + scale_local * z[1] ;
  for (i in 1:(p - 1)){
    eta[i] = z[i+1] * lambda_tilde[i] * tau;
    theta[i + 1] = eta[i] + theta[i];
  }
}

model {
  // Normal likelihood
  for (i in 1:p){
    pi_scaled[i] ~ normal(theta[i] / sigma_pi[i], 1);
  }
  
  // Priors for hyperparameters
  z ~ normal(0, 1);
  mu ~ normal(0, 1);
  
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5, 0.5);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5, 0.5);
  caux ~ inv_gamma(0.5* slab_df , 0.5* slab_df);
}
