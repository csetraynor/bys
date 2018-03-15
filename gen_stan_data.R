# !diagnostics off
##--- Gen Stan Data Null Model---#
gen_stan_data <- function(data) {

  Jobs = data$cohort

  stan_data <- list(
    Nobs = nrow(data),
    yobs = data$dfs_time,
    J = n_distinct(data$cohort),
    v = as.numeric(data$dfs_status),
    Jobs = as.numeric(Jobs)
  )
}
into_data <- gen_stan_data(md)
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))

gen_inits <- function(J) {
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1),
    
    zeta = rnorm(1),
    b = rnorm(J),
    kappa = abs(rcauchy(1, 0, 2))
  )
}
inits <- gen_inits(J = 6)
rstan::stan_rdump(ls(inits), file = "checking.init.R",
                  envir = list2env(inits))

##--- Gen Stan Data Clinical Model---#
gen_stan_data <- function(data, formula = as.formula(~1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)

  Zobs <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Zobs)
  
  Jobs = data$cohort
  
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Zobs))
      Zobs <- array(Zobs[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Zobs)
  }
  prop <- apply(Z, 2, sum) / apply(Z, 2, length) 
  Z <- sweep(Z, 2, prop) #centering
  
  stan_data <- list(
    Nobs = nrow(data),
    yobs = data$dfs_time,
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data), M)),
    v = as.numeric(data$dfs_status),
    J = n_distinct(data$cohort),
    Jobs = as.numeric(Jobs)
  )
}
Zobs <- md %>%
  model.matrix(~ stage +er +pr+ her2+menopause, data = .)
attr(Zobs, "dimnames")[[2]][-1]
into_data <- gen_stan_data(md, formula = "~ stage + er + pr + her2")
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))


gen_inits <- function(J, M) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      beta_b_raw = rnorm(M),
      
      zeta = rnorm(1),
      b = rnorm(J),
      kappa = abs(rcauchy(1, 0, 2))
    )
}
inits <- gen_inits(J = 6)
rstan::stan_rdump(ls(inits), file = "checking.init.R",
                  envir = list2env(inits))

