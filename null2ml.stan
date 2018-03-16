/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> J1; //cohort
  int<lower=0> J2; //menopause state
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int J1obs[Nobs];
  int J1cen[Ncen]; //censoring indicator
  int J2obs[Nobs];
  int J2cen[Ncen]; //censoring indicator
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  
  N = Nobs + Ncen;

  tau_al = 10.0;
  tau_mu = 10.0;
}

parameters {
  real alpha_raw;
  real mu;
  
  real zeta1;
  vector[J1] b1;
  real<lower=0> kappa1;
    
  real zeta2;
  vector[J2] b2;
  real<lower=0> kappa2;
  
}

transformed parameters {
  real<lower=0> alpha;
  vector<lower=0>[Nobs] sobs;
  vector<lower=0>[Ncen] scen;
  
  alpha = exp(tau_al * alpha_raw);
  for (n in 1:Nobs){
    sobs[n] = exp(-(mu + zeta1 + b1[J1obs[n]] + zeta1 + b1[J1obs[n]] )/alpha);
  }
  for (n in 1:Ncen){
    scen[n] = exp(-(mu + zeta1 +  b1[J1cen[n]] + zeta2 + b2[J2cen[n]]  )/alpha);
  }
}

model {
  yobs ~ weibull(alpha, sobs);
  target += weibull_lccdf(ycen | alpha, scen);

  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta1 ~ normal(0, 1);
  b1 ~ normal(0, kappa1);
  kappa1 ~ gamma_lpdf(2, .1);
  
  zeta2 ~ normal(0, 1);
  b2 ~ normal(0, kappa2);
  kappa2 ~ gamma_lpdf(2, .1);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, sobs[n]);
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, scen[n]);
  }
}