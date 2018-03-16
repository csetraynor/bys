# !diagnostics off
#Imputation#
library(mice)
md %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
#using imputation by Bayesian poly regression
md$stage[md$stage == "Unavailable"] <- NA
tmp <- as.factor(md$stage)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ pr + er +her2 ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$stage[is.na(md$stage)] <- tmp[is.na(md$stage)]
remove(tmp)

#using imputation by Bayesian poly regression
md$pr[md$pr == "Unavailable"] <- NA
tmp <- as.factor(md$pr)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ stage + er + her2 ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$pr[is.na(md$pr)] <- tmp[is.na(md$pr)]
remove(tmp)

#using imputation by Bayesian poly regression
md$er[md$er == "Unavailable"] <- NA
tmp <- as.factor(md$er)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ pr + stage +her2 ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$er[is.na(md$er)] <- tmp[is.na(md$er)]
remove(tmp)
#using imputation by Bayesian poly regression
md$her2[md$her2 == "Unavailable"] <- NA
tmp <- as.factor(md$her2)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ pr + er + stage ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$her2[is.na(md$her2)] <- tmp[is.na(md$her2)]
remove(tmp)
#using imputation by Bayesian poly regression
md$menopause[md$menopause == 4] <- NA
tmp <- as.factor(md$menopause)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ pr + er + stage + her2 ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$menopause[is.na(md$menopause)] <- tmp[is.na(md$menopause)]
remove(tmp)
