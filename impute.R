#Imputation#
md %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
md$dfs_months <- as.numeric(md$dfs_months)
md$ajcc_pathologic_tumor_stage <- as.factor(md$ajcc_pathologic_tumor_stage)


VIM::marginplot(md[c("ajcc_pathologic_tumor_stage","dfs_months")])
table(md$ajcc_pathologic_tumor_stage, useNA = "always")
#Following NPI and AJCC guidelines
md$stage <- NA
md$stage[grepl("I$|IA$|IB$",md$ajcc_pathologic_tumor_stage )] <- "1"
md$stage[grepl("II$|IIA$|IIB$",md$ajcc_pathologic_tumor_stage )] <- "2"
md$stage[grepl("IIIA$|IIIA$|IIIC$",md$ajcc_pathologic_tumor_stage )] <- "3"
md$stage[grepl("IV$|X$",md$ajcc_pathologic_tumor_stage )] <- "4"
md$stage[is.na(md$stage)] <- NA
#using imputation by Bayesian poly regression
tmp <- as.factor(md$stage)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ ajcc_nodes_pathologic_pn  + dfs_status + cancer_type_detailed +dfs_months,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$stage[is.na(md$stage)] <- tmp[is.na(md$stage)]
remove(tmp)

md$nodes <- NA
md$nodes[grepl("0",md$ajcc_nodes_pathologic_pn)] <- "1"
md$nodes[grepl("1",md$ajcc_nodes_pathologic_pn)] <- "2"
md$nodes[grepl("2|3|X",md$ajcc_nodes_pathologic_pn)] <- "3"

#--- Impute er_status and pr_status
md$erandpr <- "Negative"
md$erandpr[md$er_status_by_ihc == "Positive" & md$pr_status_by_ihc == "Positive"] <- "Positive"
md$erandpr[is.na(md$er_status_by_ihc) | md$pr_status_by_ihc == "Indeterminate"] <- NA


tmp <- as.factor(md$erandpr)
VIM::marginplot(md[c("erandpr","dfs_months")])
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~dfs_status+
                                                    cancer_type_detailed + dfs_months ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$erandpr[is.na(md$erandpr)] <- tmp[is.na(md$erandpr)]
remove(tmp)

#--- HER2 score
md$ihc_her2[is.na(md$ihc_her2) | md$ihc_her2 == "Equivocal" | md$ihc_her2 == "Indeterminate" ] <- NA
tmp <- as.factor(md$ihc_her2)

tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~dfs_status+
                                                    cancer_type_detailed + dfs_months ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$ihc_her2[is.na(md$ihc_her2)] <- tmp[is.na(md$ihc_her2)]
remove(tmp)
