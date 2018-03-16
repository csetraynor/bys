# !diagnostics off
##--- Gen Stan Data Null Model---#
gen_stan_data <- function(data) {
    J1group = data$cohort
    ind = data$dfs_recurred == 1
    
    stan_data <- list(
      Nobs = nrow(data[ind,]),
      Ncen = nrow(data[-ind,]),
      yobs = data$dfs_time[ind],
      ycen = data$dfs_time[-ind],
      J1 = n_distinct(data$cohort),
      J1obs = as.numeric(J1group[ind]),
      J1cen = as.numeric(J1group[-ind])
    )
}
# into_data <- gen_stan_data(md)
# glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits <- function(J1) {
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1),
    
    zeta1 = rnorm(1),
    b1 = rnorm(J1),
    kappa1 = abs(rcauchy(1, 0, 2))
  )
}
# inits <- gen_inits(J = 6)
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))

##--- Gen Stan Data Clinical Model---#
gen_stan_data2 <- function(data, formula = as.formula(~1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)

  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  
  Jgroup = as.integer(data$cohort)
  
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- array(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }
  prop <- apply(Z, 2, sum) / apply(Z, 2, length) 
  Z <- sweep(Z, 2, prop) #centering
  
  #Censoring indicator
  ind = data$dfs_recurred == 1
  Zobs <- Z[ind,]
  Zcen <- Z[-ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[-ind,]),
    yobs = data$dfs_time[ind],
    ycen = data$dfs_time[-ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[-ind,]), M)),
    J = n_distinct(data$cohort),
    Jobs = as.numeric(Jgroup[ind]),
    Jcen = as.numeric(Jgroup[-ind])
  )
}
# Z <- md %>%
#   model.matrix(~ stage + er + pr+ her2 + menopause , data = .)
# attr(Z, "dimnames")[[2]][-1]
into_data <- gen_stan_data2(md, formula = "~ stage + er + pr+ her2 ")
glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits2 <- function(J, M) {
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
# inits <- gen_inits(J = 6)
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))

##--- Gen Stan Data Genomic Model---#
gen_stan_data3 <- function(data, eset = NA) {
  
  Z <- t(exprs(eset))
  M <- ncol(Z)
  
  Jgroup = as.integer(data$cohort)
  
  #Censoring indicator
  ind = data$dfs_recurred == 1
  Zobs <- Z[ind,]
  Zcen <- Z[-ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[-ind,]),
    yobs = data$dfs_time[ind],
    ycen = data$dfs_time[-ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[-ind,]), M)),
    J = n_distinct(data$cohort),
    Jobs = as.numeric(Jgroup[ind]),
    Jcen = as.numeric(Jgroup[-ind])
  )
}
# Z <- md %>%
#   model.matrix(~ stage + er + pr+ her2 + menopause , data = .)
# attr(Z, "dimnames")[[2]][-1]
into_data <- gen_stan_data3(md, eset = brcaES)
glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits3 <- function(J, M) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s1_g_raw = 0.1*abs(rnorm(1)),
      tau_s2_g_raw = 0.1*abs(rnorm(1)),
      beta_g_raw = abs(rnorm(M)),
      tau1_g_raw = abs(rnorm(M)),
      tau2_g_raw = abs(rnorm(M)),
      
      zeta = rnorm(1),
      b = rnorm(J),
      kappa = abs(rcauchy(1, 0, 2))
    )
}

##--- Gen Stan Data Clinico Genomic Model---#
gen_stan_data4 <- function(data, eset = NA , formula = as.formula(~1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  
  Jgroup = as.integer(data$cohort)
  
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- array(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }
  prop <- apply(Z, 2, sum) / apply(Z, 2, length) 
  Z <- sweep(Z, 2, prop) #centering
  
  Z_g <- t(exprs(eset))
  M_g <- ncol(Z_g)
  
  Jgroup = as.integer(data$cohort)
  
  #Censoring indicator
  ind = data$dfs_recurred == 1
  Zobs <- Z[ind,]
  Zcen <- Z[-ind,]
  
  Zobs_g <- Z_g[ind,]
  Zcen_g <- Z_g[-ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[-ind,]),
    yobs = data$dfs_time[ind],
    ycen = data$dfs_time[-ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[-ind,]), M)),
    M_g = M_g,
    Zobs_g = array(Zobs_g, dim = c(nrow(data[ind,]), M_g)),
    Zcen_g = array(Zcen_g, dim = c(nrow(data[-ind,]), M_g)),
    J = n_distinct(data$cohort),
    Jobs = as.numeric(Jgroup[ind]),
    Jcen = as.numeric(Jgroup[-ind])
  )
}
# Z <- md %>%
#   model.matrix(~ stage + er + pr+ her2 + menopause , data = .)
# attr(Z, "dimnames")[[2]][-1]
# into_data <- gen_stan_data4(md, eset = brcaES, formula =  "~ stage + er + pr+ her2 ")
# glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits4 <- function(J, M, M_g) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      beta_b_raw = rnorm(M),
      
      tau_s1_g_raw = 0.1*abs(rnorm(1)),
      tau_s2_g_raw = 0.1*abs(rnorm(1)),
      
      beta_g_raw = abs(rnorm(M_g)),
      
      tau1_g_raw = abs(rnorm(M_g)),
      tau2_g_raw = abs(rnorm(M_g)),
      
      zeta = rnorm(1),
      b = rnorm(J),
      kappa = abs(rcauchy(1, 0, 2))
    )
}
rm(into_data)

