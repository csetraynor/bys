/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
data {
  int<lower=0> Nobs;
  int<lower=0> J; //cohort
  vector[Nobs] yobs;
  int Jobs[Nobs];
  int v[Nobs]; //censoring indicator
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
  
}

transformed parameters {
  real<lower=0> alpha;
  vector<lower=0>[Nobs] sobs;
  
  alpha = exp(tau_al * alpha_raw);
  for (n in 1:Nobs){
    sobs[n] = exp(-(zeta + mu + b[Jobs[n]])/alpha);
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

  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta ~ normal(0, 1);
  b ~ normal(0, kappa);
  kappa ~ gamma_lpdf(2, .1);
}
