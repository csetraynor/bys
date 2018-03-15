/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }

    return res;
  }

  vector g_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> J; //cohort
  vector[Nobs] yobs;
  int Jobs[Nobs];
  int v[Nobs]; //censoring indicator
  int<lower=0> M;
  matrix[Nobs, M] Zobs;
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  tau_al = 10.0;
  tau_mu = 10.0;
}

parameters {
  real alpha_raw;
  real mu;
  
  real zeta;
  vector[J] b;
  real<lower=0> kappa;
  
  real<lower=0> tau_s_b_raw;
  vector<lower=0>[M] tau_b_raw;
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;
  vector<lower=0>[Nobs] sobs;
  
  
  beta_b = g_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;
  alpha = exp(tau_al * alpha_raw);
  for (n in 1:Nobs){
    sobs[n] = exp(-(zeta + mu+  Zobs * beta_b + b[Jobs[n]])/alpha);
  }
}

model {
  for (i in 1:Nobs){
    if(v[i] == 1){
        yobs[i] ~ weibull(alpha, sobs[i]); 
    }else{
       target += weibull_lccdf(yobs[i] | alpha, sobs[i]); 
    }
  }
  
  beta_b_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta ~ normal(0, 1);
  b ~ normal(0, kappa);
  kappa ~ gamma_lpdf(2, .1);
}

