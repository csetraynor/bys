##Run Stan
# !diagnostics off
library(dplyr)
library(readr)
library(survival)
library(rstan)
library(loo)
library(caret)
memory.limit(5e12)
stan_file_null <- "bys/bys/nullml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stannull <- rstan::stan(stan_file_null,
                       data = gen_stan_data(md),
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 1000,
                       init = gen_inits(J1 = 6))
liknull <- loo::extract_log_lik(stannull, parameter_name = "log_lik")
loonull <- loo(liknull)
print(loonull)
saveRDS(stannull, file = "bysfit/null.rds")
rm(list = c('stannull', 'liknull'))

#---------------------------------------
##Run Clinical Stan
stan_file_clin <- "bys/bys/clinml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanclin <- rstan::stan(stan_file_clin,
                       data = gen_stan_data2(md,
                       formula = "~ stage + er + pr+ her2 "),
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 1000,
                       init = gen_inits2(J = 6, M = 12))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

likclin <- loo::extract_log_lik(stanclin, parameter_name = "log_lik")
looclin <- loo(likclin)
print(looclin)
saveRDS(stanclin, file = "bysfit/clin.rds")
rm(list = c('stanclin', 'likclin'))
compare(loonull, looclin)

#---------------------------------------
##Run Genomic Stan
stan_file_gen <- "bys/bys/genml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stangene <- rstan::stan(stan_file_gen,
                        data = gen_stan_data3(md,
                                             es = brcaES),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits3(J = 6, M = 14666))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

likgene <- loo::extract_log_lik(stangene, parameter_name = "log_lik")
loogene <- loo(likgene)
saveRDS(stangene, file = "bysfit/gene.rds")
rm(list = c('stangene', 'likgene'))


#---------------------------------------
##Run ClinicoGenomic Stan
stan_file_clingen <- "bys/bys/clingenml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanclingene <- rstan::stan(stan_file_clingen,
                        data = gen_stan_data4(md,
                        es = brcaES, formula =  "~ stage + er + pr+ her2 "),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits4(J = 6, M = 12, M_g = 14666))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

likclingene <- loo::extract_log_lik(stanclingene, parameter_name = "log_lik")
looclingene <- loo(likclingene)
saveRDS(stanclingene, file = "bysfit/clingene.rds")
rm(list = c('stanclingene', 'likclingene'))
