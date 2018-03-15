##Run Stan
# !diagnostics off
library(dplyr)
library(readr)
library(caret)
library(mice)
library(survival)
memory.limit(5e11)
stan_file <- "bys/nullml.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 1

stanfit <- rstan::stan(stan_file,
                       data = gen_stan_data(md),
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 2000,
                       init = gen_inits(J = 6))

#---------------------------------------
##Run Clinical Stan
library(dplyr)
library(rstan)
stan_file <- "bys/clinml.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 1

stanfit <- rstan::stan(stan_file,
                       data = gen_stan_data(md,
                       formula = "~ stage + er + pr+ her2 + menopause"),
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 2000,
                       init = gen_inits(J = 6, M = 16))
