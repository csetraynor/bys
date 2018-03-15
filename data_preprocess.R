# !diagnostics off
library(dplyr)
library(readr)
library(caret)
library(mice)
library(survival)
memory.limit(5e11)
#download data
sampletcga <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_sample.txt", skip = 4)
patienttcga <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_patient.txt")
gentcga <- read_tsv("brca_data/brca_tcga_pub2015/data_expression_median.txt", col_names = TRUE)

dfmeta <- read_tsv("brca_data/brca_metabric/Complete_METABRIC_Clinical_Survival_Data.txt")
samplemeta <- read_tsv("brca_data/brca_metabric/data_clinical_sample.txt", skip = 4)
patientmeta <- read_tsv("brca_data/brca_metabric/data_clinical_patient.txt")
genmeta <- read_tsv("brca_data/brca_metabric/data_expression.txt", col_names = TRUE)

# Get an easier code to read
names(sampletcga) <- tolower(names(sampletcga))
names(patienttcga) <- tolower(names(patienttcga))
names(dfmeta) <- tolower(names(dfmeta))
names(samplemeta) <- tolower(names(samplemeta))
names(patientmeta) <- tolower(names(patientmeta))

#select vars of interest and join medical data
mdtcga <- left_join(patienttcga, sampletcga) ;rm(patienttcga, sampletcga)
mdmeta <- left_join(patientmeta, samplemeta) ;rm(patientmeta, samplemeta)

mdmeta <- mdmeta %>% 
  select(patient_id, histological_subtype,
         age_at_diagnosis, tumor_stage,
         er_ihc, pr_status, her2_status, 
         inferred_menopausal_state,
         chemotherapy,  radio_therapy,
         threegene, cohort) %>%
  rename(histology = histological_subtype,
         age = age_at_diagnosis,
         stage = tumor_stage,
         er = er_ihc,
         pr = pr_status,
         her2 = her2_status,
         menopause = inferred_menopausal_state,
         radio = radio_therapy,
         chemo = chemotherapy)

dfmeta <- dfmeta %>%
  select(patient_id, time, status, lymph_nodes_positive) %>% 
  rename(dfs_time = time,
         dfs_status = status,
         nodes = lymph_nodes_positive)
mddfmeta <- inner_join(dfmeta, mdmeta);rm(dfmeta, mdmeta)

mdtcga <- mdtcga %>%  
  select(sample_id, dfs_time, dfs_status, 
         ajcc_nodes_pathologic_pn,
         histological_diagnosis,
         age, ajcc_pathologic_tumor_stage,
         er_status_by_ihc, pr_status_by_ihc, ihc_her2, 
         menopause_status,  pharmaceutical_tx_adjuvant,
         radiation_treatment_adjuvant) %>% 
  rename(patient_id = sample_id,
         histology = histological_diagnosis,
         nodes = ajcc_nodes_pathologic_pn,
         stage = ajcc_pathologic_tumor_stage,
         er = er_status_by_ihc,
         pr = pr_status_by_ihc,
         her2 = ihc_her2,
         menopause = menopause_status,
         chemo = pharmaceutical_tx_adjuvant,
         radio = radiation_treatment_adjuvant) %>%
  mutate( threegene = NA , cohort = "6",dfs_status = (dfs_status == "Recurred/Progressed"))
assertthat::assert_that(all(colnames(mdtcga) == colnames(mddfmeta)))
md <- rbind(mddfmeta , mdtcga);rm(mddfmeta,mdtcga)
glimpse(md)
#---- Data Cleaning ----#
#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    return(x)
  } else {
    ifelse(x == " " | x == "[Not Available]" | x == "[Discrepancy]" | x == "Indeterminate"| x == "Equivocal", NA, x)
  }
}
md <- md %>% dplyr::mutate_all(funs(convert_blank_to_na)) %>%
  mutate_at(vars("age"), funs(as.numeric))

md %>% VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
            sortVars = TRUE, sortCombs = TRUE, plot = TRUE, only.miss = FALSE)
md %>% 
  filter(is.na(dfs_status) | dfs_status == "" |dfs_time <= 0 | is.na(dfs_time)) %>%
  select(dfs_time, dfs_status) %>%
  glimpse
#--- Remove Missing Obs ---#
clinical_data <- md
md <- md %>% 
  filter(!is.na(dfs_status) , dfs_status != "" , dfs_time > 0 , !is.na(dfs_time))
#confirm 78 fewer observations
assertthat::assert_that(nrow(md) == (nrow(clinical_data) - 90))
rm(clinical_data)
#--- Distribution event times ---#
md <- md %>%
  arrange(dfs_time)
md <- md %>% mutate(dfs_status = dfs_status == "1")
md %>%
  ggplot(aes(x = dfs_time,
             group = dfs_status,
             colour = dfs_status,
             fill = dfs_status)) +
  geom_density(alpha = 0.5)
mle.surv <- survfit(Surv(dfs_time, dfs_status) ~ cohort,
                    data = md )
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM estimate for BCs Cohort')

#---- Pipe ---#

glimpse(md)
md$er[md$er == "Positive" |md$er == "+" ] <- "pos"
md$er[md$er == "Negative" | md$er == "-"] <- "neg"
md$er[md$er == "Indeterminate" | is.na(md$er)] <- "Unavailable" 

md$pr[md$pr == "Positive" |md$pr == "+"] <- "pos"
md$pr[md$pr == "Negative" | md$pr == "-"] <- "neg"
md$pr[md$pr == "Indeterminate" | is.na(md$pr)] <- "Unavailable" 

md$her2[md$her2 == "Positive" |md$her2 == "+" ] <- "pos"
md$her2[md$her2 == "Negative" | md$her2 == "-"] <- "neg"
md$her2[md$her2 == "Indeterminate" | md$her2 == "Equivocal"  | is.na(md$her2) ] <-"Unavailable"  

md$menopause[md$menopause == "Peri (6-12 months since last menstrual period)" |
               md$menopause == "Indeterminate (neither Pre or Postmenopausal)"] <- "xperi"
md$menopause[md$menopause == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"] <- "post"
md$menopause[md$menopause == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"] <- "pre"
md$menopause[is.na(md$menopause)] <- "Unavailable"

md$stage[grepl("0$",md$stage )] <- "0"
md$stage[grepl("I$|IA$|IB$",md$stage )] <- "1"
md$stage[grepl("II$|IIA$|IIB$",md$stage )] <- "2"
md$stage[grepl("IIIA$|IIIA$|IIIC$",md$stage )] <- "3"
md$stage[grepl("IV$",md$stage )] <- "4"
md$stage[grepl("X$", md$stage)] <- "null"
md$stage[is.na(md$stage)] <- "Unavailable"


#--- Gene matrix preprocess ----- #
source("https://bioconductor.org/biocLite.R")
library(Biobase)
gendat <- inner_join(genmeta, gentcga, by = "Hugo_Symbol");rm(genmeta, gentcga)
glimpse(gendat)
gene_names <- gendat %>% select(Hugo_Symbol)  #grap gene names
gendat <- gendat %>% select(intersect(colnames(gendat), md$patient_id))# get intersection btw clinical and gene values
sample_names <- colnames(gendat) # get sample names
sum(is.na(gendat))
#Convert to expression set
md <- as.data.frame(md %>% filter(patient_id %in% sample_names) %>% slice(match(sample_names, patient_id))) ; row.names(md) <- md$patient_id#requires classical data frame
x <-  as.matrix(gendat) ;colnames(x) <- rownames(md); 
gene_names <- as.data.frame(gene_names);  rownames(gene_names) <- gene_names %>% unlist
brcaES <- Biobase::ExpressionSet(x,
                                  phenoData = as(md, "AnnotatedDataFrame"),
                                 featureData = as(gene_names, "AnnotatedDataFrame"))
assertthat::assert_that(all(md$patient_id == brcaES$patient_id))
rm(list = "x")
gene_names <- gene_names %>% unlist
#Imputation using MSnbase
require(MSnbase)
brcaMSN <- MSnbase::as.MSnSet.ExpressionSet(brcaES)
brcaMSN <- MSnbase::impute(brcaMSN, method = "knn")
Biobase::exprs(brcaES) <- MSnbase::exprs(brcaMSN)
rm(brcaMSN)
sum(is.na(Biobase::exprs(brcaES)))

##Perform empirical Bayes to find differential gene expressions
library(limma)
fit <- limma::eBayes(limma::lmFit(brcaES))
volcanoplot(fit)
# toptable(fit)
# # 
rm(fit)
#Center and scale 
preProcgm <-  caret::preProcess(t(exprs(brcaES)), method = c("center", "scale")) 
brcaES <- predict(preProcgm, t(exprs(brcaES))) 
rm(preProcgm)

# # ### x is the input data. This function replaces the top 'perc' percent
# # ### with the value 'rp'. 
# # 
# subset.top = function(x, perc, eset)
# {
#   qnt = quantile(x$lods,1-perc)
#   w = which(fit$lods >= qnt)
#   return(eset[w])
# }
# brcaES = subset.top(x = fit, perc = .2 , eset = brcaES)
# volcanoplot(fit)
# fit_gene_names = rownames(brcaEStest)
# rm(f
