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

  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
  }
}
data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> J; //cohort
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int Jobs[Nobs];
  int Jcen[Ncen]; //censoring indicator
  int<lower=0> M;
  matrix[Nobs, M] Zobs;
  matrix[Ncen, M] Zcen;
  int<lower=0> M_g;
  matrix[Nobs, M_g] Zobs_g;
  matrix[Ncen, M_g] Zcen_g;
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  int<lower=0> nu;
  tau_al = 10.0;
  tau_mu = 10.0;
  N = Nobs + Ncen;
  nu = 3;
}

parameters {
  real alpha_raw;
  real mu;
  
  real zeta;
  vector[J] b;
  real<lower=0> kappa;
  
  real<lower=0> tau_s1_g_raw;
  real<lower=0> tau_s2_g_raw;
  vector<lower=0>[M_g] tau1_g_raw;
  vector<lower=0>[M_g] tau2_g_raw;
  vector[M_g] beta_g_raw;
  
  real<lower=0> tau_s_b_raw;
  vector<lower=0>[M] tau_b_raw;
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M_g] beta_g;
  vector<lower=0>[Nobs] sobs;
  vector<lower=0>[Ncen] scen;
  vector[M] beta_b;
  
  
  beta_g = hs_prior_lp(tau_s1_g_raw, tau_s2_g_raw, tau1_g_raw, tau2_g_raw, nu) .* beta_g_raw;
  beta_b = g_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;
  alpha = exp(tau_al * alpha_raw);
  
  for (n in 1:Nobs){
    sobs[n] = exp(-( mu + zeta  + Zobs[n,] * beta_b + Zobs_g[n,] * beta_g + b[Jobs[n]])/alpha);
  }
  for (n in 1:Ncen){
    scen[n] = exp(-( mu + zeta  + Zcen[n,] * beta_b + Zcen_g[n,] * beta_g + b[Jcen[n]])/alpha);
  }
}

model {
  yobs ~ weibull(alpha, sobs);
  target += weibull_lccdf(ycen | alpha, scen);
  
  beta_g_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta ~ normal(0, 1);
  b ~ normal(0, kappa);
  kappa ~ gamma_lpdf(2, .1);
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


