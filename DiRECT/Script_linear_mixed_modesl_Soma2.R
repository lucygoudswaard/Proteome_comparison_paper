#### Running linear mixed model script on soma2
##### Script to run analyses for manuscript
## load packages
library(openxlsx)
library(zip)
library(data.table)
library(systemfit)
library(ggplot2)
library(devtools)
library(dplyr)
library(cowplot)
library(reshape2)
library(ggbiplot)
library(tidyr)
library(purrr)
library(tidyverse)
library(tidyr)
library(broom)
library(generics)
library(glue)
library(utiles)
library(mice)
library(zoo)
library(psych)
library(BiocManager)
#library(doppelgangR)
library(readr)
library(MetaboQC)
library(lme4)
library(readxl)

## Read in phenotype data
all_data <- readRDS(file = "P:/DiRECT_fu1_20171113.rds")

#names(all_data)
print(paste("Number of individuals in trial dataset:", nrow(all_data)))

# retain variables for model
model_data <- all_data[,c("id","site","centre","treat","list.size","age","sex","bmi.b","bmi","bmi.change","weight.b","weight","weight.change", "sbp.b",
                          "dbp.b", "diabdur", "n.anti.diab.b", "hba1c.mmol.mol.b")]
#dim(model_data)
#head(model_data)
#tail(model_data)

# rename weight and bmi endpoint var so .e added to signify endpoint
names(model_data)[9] <- "bmi.e"
names(model_data)[12] <- "weight.e"

# check class of factors for model
class(model_data$centre)
class(model_data$treat)
class(model_data$list.size)
class(model_data$site)
class(model_data$age)
class(model_data$sex)

# reformat as required
model_data$site <- as.factor(model_data$site)

# retain clinical bloods variables for analysis
bloods_data <- all_data[,c(1,126,127,129,130,139,140,141,152,153,154)]

print("Blood variables retained for analysis alongside metabolites:")
names(bloods_data)

# rename var - remove '.'
names(bloods_data)[2] <- "glucose_mmol_l.b"
names(bloods_data)[3] <- "glucose_mmol_l.e"
names(bloods_data)[4] <- "insulin_uu_ml.b"
names(bloods_data)[5] <- "insulin_uu_ml.e"
names(bloods_data)[6] <- "hdl_mmol_l.b"
names(bloods_data)[9] <- "hdl_mmol_l.e"
names(bloods_data)[7] <- "trig_mmol_l.b"
names(bloods_data)[10] <- "trig_mmol_l.e"
names(bloods_data)[8] <-"chol_mmol_l.b"
names(bloods_data)[11] <- "chol_mmol_l.e"

## Merge bloods_data and model_data
basic_phenos <- merge(model_data, bloods_data, by = "id")

## create dataframe to do T2D clustering
## subset other variables required for clustering - need zscore age, zscore bmi, hba1c zscore, HOMA2-B szore, HOMA-IR zscore
## write excel file with ids, insulin, glucose, age, bmi, hba1c that can be entered into the calculator
homa_calc_data <- all_data[,c('id', 'insulin.uu.ml.b', 'glucose.mmol.l.b', 'hba1c.mmol.mol.b', 'age', 'bmi.b' )]
#write.table(homa_calc_data, file = 'U:/Results/Clustering_data.txt', sep = '\t', col.names = T, row.names = F)

all_data$HOMAIR <- (all_data$insulin.uu.ml.b * all_data$glucose.mmol.l.b) / 405
all_data$HOMAB <- (20* all_data$insulin.uu.ml.b) / (all_data$glucose.mmol.l.b - 3.5)
w <- which(all_data$HOMAB > 4422)
all_data$HOMAB[w] <- NA
cluster_dataframe <- all_data[, c('age', 'bmi.b', 'hba1c.mmol.mol.b', 'HOMAIR' , 'HOMAB')]




## Column ID match
ColumnInfo <- fread("P:/19022021/ColumnInfo.csv", sep = ",", check.names = F, header = T)
ColumnInfo <- as.data.frame(ColumnInfo)

## Remove mouse proteins from protein info dataframe
w <- which(ColumnInfo[8,] == "Mouse")
ColumnInfo <- ColumnInfo[,-w]

## Dataframe matching proteins to merge to regression results
proteinmatch <- t(ColumnInfo)
colnames(proteinmatch) <- proteinmatch[1,]
proteinmatch <- proteinmatch[-1,]
proteinmatch <- as.data.frame(proteinmatch)
proteinmatch <- proteinmatch[, c("SomaId", "TargetFullName", "Target", "UniProt", "EntrezGeneID", "EntrezGeneSymbol")]
proteinmatch$UniProt <- as.character(proteinmatch$UniProt)
idsoma <- rownames(proteinmatch)
proteinmatch <- cbind(idsoma, proteinmatch)
rownames(proteinmatch) <- c()

## Read in soma2 dataset which has the extra normalization step
soma2 <- readRDS(file = "P:/Somalogic/soma2.rds")
##Remove those that say visit 0 (should just be visit 1 and 2)
w <- which(soma2$visit == "v0")
soma2 <- soma2[-w,]
# Order dataframe by ID
soma2 <- soma2[order(soma2$id)]
soma2$visit <- as.factor(soma2$visit)
soma2$PlateId <- as.factor(soma2$PlateId)
phenos <- soma2[,c(4635, 4636, 1:33)]
phenos$idlong <- paste0(phenos$id, '_', phenos$visit)
## Drop all "SeqId." from cols 34:4634
colnames(soma2)[34:4634] <- sub("SeqId.", "", colnames(soma2)[34:4634])
length(unique(soma2$id)) ## 302 individuals

metaboprep_soma2 <- soma2
metaboprep_soma2 <- metaboprep_soma2[,c(4635, 4636, 34:4634)]
metaboprep_soma2$id_long <- paste0(as.character(metaboprep_soma2$id), '_', metaboprep_soma2$visit)
metaboprep_soma2 <- metaboprep_soma2[,c(4604, 3:4603)]

## Write file for metaboprep 
#write.table(metaboprep_soma2, file = "U:/metaboprep/metaboqc_data_soma2.txt", sep = "\t", col.names = T, row.names = F)

## Reread in metaboprep data and replace soma data 
soma2_qc <- read.table(file = 'U:/metaboprep/metaboprep_release_2022_11_30/filtered_data/DiRECT_SomaLogic_2_2022_11_30_Filtered_metabolite_data.txt', sep = '\t', header = T)
idlong <- rownames(soma2_qc)
soma2_qc <- cbind(idlong, soma2_qc)
colnames(soma2_qc)[2:4602] <- colnames(soma2)[34:4634]
soma2_qc$idlong <- as.character(soma2_qc$id)


## Merge qc soma2 back to 
full_qc <- merge(phenos, soma2_qc, by = "idlong", all.y = T)
full_qc <- full_qc[,-1]
length(unique(full_qc$id)) ## 300 individuals

## Remove .1 after Id names
full_qc$id <- gsub(pattern = "\\..*", replacement = "", full_qc$id)

## Natural log transformation
# full_qc <- as.data.frame(full_qc)
# for (i in 36:4634){
#   full_qc[,i] <- log2(full_qc[,i])
# }

## Merge in age and sex
agesex <- basic_phenos[,c("id", "age", "sex")]
full_qc <- merge(agesex, full_qc, by = "id")
length(unique(full_qc$id)) ## 292

## rnt function
rnt <- function(x){
  qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))
}

## Change soma1 to wide format
soma1.wide <- pivot_wider(
  full_qc,
  names_from = visit,
  names_sep = "_",
  values_from = 5:4638
)

## change _v1 to .b and _v2 to .e
colnames(soma1.wide)[69:9270] <- gsub("_v1", ".b", colnames(soma1.wide)[69:9270])
colnames(soma1.wide)[69:9270] <- gsub("_v2", ".e", colnames(soma1.wide)[69:9270])

## combine basic_phenos, soma1.wide (ca) 
#basic_phenos <- basic_phenos[, -c(6:7)]
## This dataset can be used if using protein variables split into timepoints
data <- merge(basic_phenos, soma1.wide, by = "id", all.x = FALSE)
data$id <- as.character(data$id)

## Summaries of baseline data
colnames <- colnames(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                                "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                                "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")])
means <- apply(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                          "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                          "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")], 2, function(x) { mean(na.omit(x))})
SDs <- apply(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                        "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                        "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")], 2, function(x) { sd(na.omit(x))})
n <- apply(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                      "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                      "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")], 2, function(x) { length(na.omit(x))})
minimum <- apply(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                            "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                            "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")], 2, function(x) { min(na.omit(x))})
maximum <- apply(data[, c("age.x", "bmi.b", "bmi.e", "bmi.change", "glucose_mmol_l.b", "glucose_mmol_l.e",
                            "insulin_uu_ml.b", "insulin_uu_ml.e", "hdl_mmol_l.b", "hdl_mmol_l.e", "trig_mmol_l.b",
                            "trig_mmol_l.e", "chol_mmol_l.b", "chol_mmol_l.e")], 2, function(x) { max(na.omit(x))})
sumstats <- cbind(colnames, means, SDs, n, minimum, maximum)

datsum <- as.data.frame(sumstats)
datsum

## Use this if only wanting complete case data
#data <- data[complete.cases(data[,c(90:9291)]),]

## Use "data" dataframe to create baseline characteristics table
data$treat <- as.factor(data$treat)
data$centre <- as.factor(data$centre)
data$list.size <- as.factor(data$list.size)
data$site <- as.factor(data$site)

## Fishers exact or chi squared test to see if treatment is associated with sex, centre, list.size, site....
sex <- chisq.test(data$sex.x, data$treat)
sex ## P=0.2844
table(data$sex.x, data$treat) ## Control F 55, M = 91, Treat F 65, M = 81
centre <- chisq.test(data$centre, data$treat)
centre ## p = 0.0001788
table(data$centre, data$treat) ## Control S 114, T 32, Treat S 83. T 63
age <- t.test(data$age.x ~ data$treat) ## Control ##56.2, Treat 53.0
age ## 0.0002329
length(na.omit(data$age.x))
control <- which(data$treat == "Control")
sdagecontrol <- sd(data$age.x[control]) ## 7.06
treat <- which(data$treat == "Intervention") 
sdagetreat <- sd(data$age.x[treat]) ## 7.53
listsize <- chisq.test(data$list.size, data$treat)
listsize ## p=0.05936
table(data$list.size, data$treat) ## Control > 5700 73, <5700 73, Int > 5700 90, <5700 56
site <- chisq.test(data$site, data$treat)
site$p.value ## 4.30e-37
basebmi <- t.test(data$bmi.b ~ data$treat)
length(na.omit(data$bmi.b))
basebmi ## p = 0.09877
sdbmicontrol <- sd(data$bmi.b[control]) ## 4.27
sdbmitreat <- sd(data$bmi.b[treat]) ## 4.56
basesbp <- t.test(data$sbp.b ~ data$treat)
basesbp ## 0.02595
sdsbpcontrol <- sd(data$sbp.b[control]) ## 15.8
sdsbptreat <- sd(data$sbp.b[treat]) ## 17.38
length(na.omit(data$sbp.b))
bmi_change <- t.test(data$bmi.change ~ data$treat)
bmi_change ## mean control = -0.34, mean intervention -3.50, p=3.36x10-25
sdbmichangecontrol <- sd(na.omit(data$bmi.change[control])) ## 1.32
sdbmichangetreat <- sd(na.omit(data$bmi.change[treat])) ## 2.77

## Repeat comparisons for sex
centre <- chisq.test(data$centre, data$sex.x)
centre ## p = 0.7111
table(data$centre, data$sex.x) 

treattest <- chisq.test(data$treat, data$sex.x)
treattest ## p = 0.28
table(data$treat, data$sex.x)


age <- t.test(data$age.x ~ data$sex.x)
age ## 0.95
Female <- which(data$sex.x == "Female")
sdageFemale <- sd(data$age.x[Female]) 
Male <- which(data$sex.x == "Male")
sdageMale <- sd(data$age.x[Male]) 
listsize <- chisq.test(data$list.size, data$sex.x)
listsize ## p=0.90
table(data$list.size, data$sex.x)
basebmi <- t.test(data$bmi.b ~ data$sex.x)
basebmi ## p = 0.12
sdbmiFemale <- sd(data$bmi.b[Female]) 
sdbmiMale <- sd(data$bmi.b[Male]) 
basesbp <- t.test(data$sbp.b ~ data$sex.x)
basesbp ## 0.10
sdsbpFemale <- sd(data$sbp.b[Female]) 
sdsbpMale <- sd(data$sbp.b[Male]) 

ndiabmed <- t.test(data$n.anti.diab.b ~ data$sex.x)
ndiabmed ## 0.01805
table(data$n.anti.diab.b, data$sex.x)

basehba1c <- t.test(data$hba1c.mmol.mol.b ~ data$sex.x)
basehba1c ## 0.30
sdhba1cFemale <- sd(data$hba1c.mmol.mol.b[Female]) ## 15.9
sdhba1cMale <- sd(data$hba1c.mmol.mol.b[Male]) ## 9.6

baseglucose <- t.test(data$glucose_mmol_l.b ~ data$sex.x)
baseglucose ## 0.63
sdglucoseFemale <- sd(na.omit(data$glucose_mmol_l.b[Female]))## 3.4
sdglucoseMale <- sd(na.omit(data$glucose_mmol_l.b[Male])) ## 2.96

baseinsulin <- t.test(data$insulin_uu_ml.b ~ data$sex.x)
baseinsulin ## 0.61
sdinsulinFemale <- sd(na.omit(data$insulin_uu_ml.b[Female])) ## 18.8
sdinsulinMale <- sd(na.omit(data$insulin_uu_ml.b[Male])) ## 12.4

basechol <- t.test(data$chol_mmol_l.b ~ data$sex.x)
basechol ## 0.0014
sdcholFemale <- sd(na.omit(data$chol_mmol_l.b[Female])) #1.20
sdcholMale <- sd(na.omit(data$chol_mmol_l.b[Male])) ## 0.93

basehdl <- t.test(data$hdl_mmol_l.b ~ data$sex.x)
basehdl ## 4.1e-5
sdhdlFemale <- sd(na.omit(data$hdl_mmol_l.b[Female]))# 0.24
sdhdlMale <- sd(na.omit(data$hdl_mmol_l.b[Male])) ## 0.26

basetrig <- t.test(data$trig_mmol_l.b ~ data$sex.x)
basetrig # 0.07
sdtrigFemale <- sd(na.omit(data$trig_mmol_l.b[Female])) ## 0.84
sdtrigMale <- sd(na.omit(data$trig_mmol_l.b[Male])) ## 1.74

basediabdur <- t.test(data$diabdur ~ data$sex.x)
basediabdur ## 0.84
sddiabFemale <- sd(data$diabdur[Female]) ## 1.56
sddiabMale <- sd(na.omit(data$diabdur[Male])) ## 1.61

ndiabmed <- t.test(data$n.anti.diab.b ~ data$sex.x)
ndiabmed 
sdndiabmedFemale <- sd(na.omit(data$n.anti.diab.b[Female]))
sdndiabmedMale <- sd(na.omit(data$n.anti.diab.b[Male]))

listsize <- chisq.test(data$list.size, data$sex.x)
listsize ## p=0.90
table(data$list.size, data$sex.x)


# Max and min for continuous data
data <- as.data.frame(data)
summary_data <- data[, c('bmi.b', 'sbp.b', "n.anti.diab.b", 'glucose_mmol_l.b', 'insulin_uu_ml.b', 'chol_mmol_l.b', 'hdl_mmol_l.b', 
                         'trig_mmol_l.b', 'diabdur') ]
summary_data$n.anti.diab.b <- as.numeric(summary_data$n.anti.diab.b)
cnames <- colnames(summary_data)
for (i in 1:9) {
  var <- cnames
  mean[i] <- mean(summary_data[,i], na.rm = T)
  min[i] <- min(summary_data[,i], na.rm = T)
  max[i] <- max(summary_data[,i], na.rm = T)
  summaries <- as.data.frame(cbind(var, mean, min, max))
}

ndiabmed <- t.test(data$n.anti.diab.b ~ data$treat)
ndiabmed ## 0.5978
table(data$n.anti.diab.b, data$treat)
## mean control = 1.09, mean intervention = 1.14
length(na.omit(data$n.anti.diab.b[control])) ##146
length(na.omit(data$n.anti.diab.b[treat])) ##146
sd(na.omit(data$n.anti.diab.b[control]))
sd(na.omit(data$n.anti.diab.b[treat]))
   
basehba1c <- t.test(data$hba1c.mmol.mol.b ~ data$treat)
basehba1c ## 0.13
sdhba1ccontrol <- sd(data$hba1c.mmol.mol.b[control]) ## mean 58.1, SD 11.6
length(na.omit(data$hba1c.mmol.mol.b[control])) ## 146
sdhba1ctreat <- sd(data$hba1c.mmol.mol.b[treat]) ## mean 60.4, SD 13.8
length(na.omit(data$hba1c.mmol.mol.b[treat])) ##146

baseglucose <- t.test(data$glucose_mmol_l.b ~ data$treat)
baseglucose ## 0.2332
sdglucosecontrol <- sd(data$glucose_mmol_l.b[control]) ## mean 8.80, SD 2.55
sdglucosetreat <- sd(na.omit(data$glucose_mmol_l.b[treat])) ## mean 9.22 SD 3.29
length(na.omit(data$glucose_mmol_l.b[control])) ## 146
length(na.omit(data$glucose_mmol_l.b[treat])) ## 144


baseinsulin <- t.test(data$insulin_uu_ml.b ~ data$treat)
baseinsulin ## 0.2597
sdinsulincontrol <- sd(data$insulin_uu_ml.b[control]) #mean 22.62 SD 13.75
sdinsulintreat <- sd(na.omit(data$insulin_uu_ml.b[treat])) # mean 24.53 SD 14.97
length(na.omit(data$insulin_uu_ml.b[control])) ## 146
length(na.omit(data$insulin_uu_ml.b[treat])) ## 144


basechol <- t.test(data$chol_mmol_l.b ~ data$treat)
basechol ## 0.9086
sdcholcontrol <- sd(data$chol_mmol_l.b[control]) ## mean 4.32 SD 1.23
sdcholtreat <- sd(na.omit(data$chol_mmol_l.b[treat])) ## mean 4.34 SD 1.15
length(na.omit(data$chol_mmol_l.b[control])) ## 146
length(na.omit(data$chol_mmol_l.b[treat])) ## 143


basehdl <- t.test(data$hdl_mmol_l.b ~ data$treat)
basehdl ## 0.009499
sdhdlcontrol <- sd(data$hdl_mmol_l.b[control]) ## mean 1.16 SD 0.31
sdhdltreat <- sd(na.omit(data$hdl_mmol_l.b[treat])) ## mean 1.08 SD 0.25
length(na.omit(data$hdl_mmol_l.b[control])) ## 146
length(na.omit(data$hdl_mmol_l.b[treat])) ## 143


basetrig <- t.test(data$trig_mmol_l.b ~ data$treat)
basetrig ## 0.3375
sdtrigcontrol <- sd(data$trig_mmol_l.b[control]) ## mean 1.94 SD 0.93
sdtrigtreat <- sd(na.omit(data$trig_mmol_l.b[treat])) ## mean 2.07 SD 1.37
length(na.omit(data$trig_mmol_l.b[control])) ## 146
length(na.omit(data$trig_mmol_l.b[treat])) ## 143


basediabdur <- t.test(data$diabdur ~ data$treat)
basediabdur ## 0.8147
sddiabcontrol <- sd(data$diabdur[control]) ## mean 2.99 SD 1.75
sddiabtreat <- sd(na.omit(data$diabdur[treat])) ## mean 3.04 SD 1.68
length(na.omit(data$diabdur[control])) ## 146
length(na.omit(data$diabdur[treat])) ## 146


## W statistic of BMI
lshap <- lapply(data[,c(8:28)], shapiro.test)
lres <- sapply(lshap, `[`, "statistic")
w <- unlist(lres)
hist(data$bmi.b) ## 0.968
hist(data$bmi.e) ## 0.978


## rnt protein data (use cols 38:4638 if using the log2 and age/sex adjustment steps)
full_qc <- as.data.frame(full_qc)
for (i in 38:4638) {
  full_qc[,i] <- rnt( full_qc[,i]) }

## Adjust variables for age, sex
for (i in 38:4638){
  results <- lm(full_qc[,i] ~ full_qc$age + full_qc$sex )
  full_qc[,i] <- resid(results)
}

length(unique(full_qc$id))

## Check distributions of transformed proteomic data 
lshap <- lapply(full_qc[, 38:4638], shapiro.test)
lres <- sapply(lshap, `[`, "statistic")
w <- unlist(lres)
W_statistic_by_timepoint <- w
hist(W_statistic_by_timepoint)

## While in long format: lmer for protein ~ time
treatment <- basic_phenos[,c("id", "treat", "centre", "list.size")]
glm_data <- merge(x= full_qc, y=treatment, by = "id", all.x=TRUE) 

## How many participants are in this analysis?
length(unique(glm_data$id)) ## 292

## Run the steps below for complete case linear mixed model (this created the saved excel file lmer_protein_time_30.9.22_cc,xlsx)
# x <- table(glm_data$id)
# x <- as.data.frame(x)
# w <- which(x[,2] == 1)
# not_repeated <- x[,1][w]
# w <- which(glm_data$id %in% not_repeated)
# glm_data <- glm_data[-w,]

prot_time_results_soma2 <- as.data.frame(matrix(ncol=14, nrow=4602))
colnames(prot_time_results_soma2) <- c("Protein", "Subject_variance", "Subject_SD", "Intervention_beta", "Intervention_SE", "Intervention_t", "V2_beta",
                                 "V2_SE", "V2_t", "Intervention*v2_Beta", "Intervention*V2_SE", "Intervention*V2_t", "N", "P")
for (i in 38:4638){
  pt <- lmer( glm_data[,i] ~ glm_data$treat * glm_data$visit + glm_data$centre + glm_data$list.size + (1 | glm_data$id))
  protein_time <- summary(lmer( glm_data[,i] ~ glm_data$treat * glm_data$visit + + glm_data$centre + glm_data$list.size + (1 | glm_data$id)))
  Stdev <- as.data.frame(VarCorr(pt))
  FitA <- lmer( glm_data[,i] ~ glm_data$treat * glm_data$visit + glm_data$centre + glm_data$list.size + (1 | glm_data$id), REML = FALSE)
  # Fit null is Fit B
  FitB <-  lmer( glm_data[,i] ~ glm_data$treat + glm_data$visit + glm_data$centre + glm_data$list.size + (1 | glm_data$id), REML = FALSE)
  compare <- anova(FitA, FitB)
  prot_time_results_soma2[i-36,1] <- names(glm_data)[i]
  prot_time_results_soma2[i-36,2] <- Stdev[1,4]
  prot_time_results_soma2[i-36,3] <- Stdev[1,5]
  prot_time_results_soma2[i-36,4] <- protein_time$coefficients[2,1]
  prot_time_results_soma2[i-36,5] <- protein_time$coefficients[2,2]
  prot_time_results_soma2[i-36,6] <- protein_time$coefficients[2,3]
  prot_time_results_soma2[i-36,7] <- protein_time$coefficients[3,1]
  prot_time_results_soma2[i-36,8] <- protein_time$coefficients[3,2]
  prot_time_results_soma2[i-36,9] <- protein_time$coefficients[3,3]
  prot_time_results_soma2[i-36,10] <- protein_time$coefficients[6,1]
  prot_time_results_soma2[i-36,11] <- protein_time$coefficients[6,2]
  prot_time_results_soma2[i-36,12] <- protein_time$coefficients[6,3]
  prot_time_results_soma2[i-36,13] <- protein_time$devcomp$dims[5]
  ### This line below is wrong
  prot_time_results_soma2[i-36,14] <- compare[2,8]
  
}
prot_time_results_soma2 <- merge(prot_time_results_soma2, proteinmatch, by.x = "Protein", by.y = "idsoma", all.x = T, all.y = F, sort = F)
#write.xlsx(prot_time_results_soma2, file = "U:/Results/lmer_protein_time_soma2_30.11.22.xlsx")

# ## Stratify by linear mixed model by group - intervention arm only
# w <- which(glm_data$treat == "Intervention")
# lmer_dat_intervention <- glm_data[w,] 
# prot_time_intervention_soma2 <- as.data.frame(matrix(ncol=8, nrow=4602))
# colnames(prot_time_intervention_soma2) <- c("Protein", "Subject_variance", "Subject_SD", "V2_beta",
#                                       "V2_SE", "V2_t", "N", "P")
# for (i in 37:4637){
#   pt <- lmer( lmer_dat_intervention[,i] ~ lmer_dat_intervention$visit + lmer_dat_intervention$centre + lmer_dat_intervention$list.size +  (1 | lmer_dat_intervention$id))
#   protein_time <- summary(lmer( lmer_dat_intervention[,i] ~ lmer_dat_intervention$visit + lmer_dat_intervention$centre + lmer_dat_intervention$list.size +  (1 | lmer_dat_intervention$id)))
#   FitA <- lmer( lmer_dat_intervention[,i] ~ lmer_dat_intervention$visit + lmer_dat_intervention$centre + lmer_dat_intervention$list.size +  (1 | lmer_dat_intervention$id), REML = FALSE)
#   FitB <-  lmer( lmer_dat_intervention[,i] ~ lmer_dat_intervention$centre + lmer_dat_intervention$list.size +  (1 | lmer_dat_intervention$id), REML = FALSE)
#   compare <- anova(FitA, FitB)
#   Stdev <- as.data.frame(VarCorr(pt))
#   prot_time_intervention_soma2[i-36,1] <- names(lmer_dat_intervention)[i]
#   prot_time_intervention_soma2[i-36,2] <- Stdev[1,4]
#   prot_time_intervention_soma2[i-36,3] <- Stdev[1,5]
#   prot_time_intervention_soma2[i-36,4] <- protein_time$coefficients[2,1]
#   prot_time_intervention_soma2[i-36,5] <- protein_time$coefficients[2,2]
#   prot_time_intervention_soma2[i-36,6] <- protein_time$coefficients[2,3]
#   prot_time_intervention_soma2[i-36,7] <- 119 #protein_time$devcomp$dims[5]
#   prot_time_intervention_soma2[i-36,8] <- compare[2,8]
# }
# 
# prot_time_intervention_soma2 <- merge(prot_time_intervention_soma2, proteinmatch, by.x = "Protein", by.y = "idsoma", all.x = T, all.y = F, sort = F)
# 
# 
# full_qc$PlateId <- as.factor(full_qc$PlateId)
# regression_results5 <- as.data.frame(matrix(ncol = 7, nrow =4601))
# colnames(regression_results5) <- c("Protein", "Beta", "Pval", "SE", "AdjR", "N", "Fstat")
# for (i in 38:4638){
#   results <- summary(lm(full_qc[,i] ~ full_qc$PlateId))
#   Fstat <- as.data.frame(results$fstatistic)
#   a <- anova(lm(full_qc[,i] ~ full_qc$PlateId))
#   eta <- (a[,2]/sum(a[,2]))[1]
#   regression_results5[i-36,1] <-colnames(full_qc)[i]
#   regression_results5[i-36,2] <- results$coefficients[2,1]
#   regression_results5[i-36,3] <- results$coefficients[2,4]
#   regression_results5[i-36,4] <- results$coefficients[2,2]
#   regression_results5[i-36,5] <- results$adj.r.squared
#   regression_results5[i-36,6] <- length(residuals(results))
#   regression_results5[i-36,7] <- Fstat[1,1]
# }
# 
# regression_results5 <- regression_results5[order(regression_results5[,"Pval"]),]
# regression_results5$Protein <- gsub(pattern = '^.*\\.', replacement = "",regression_results5$Protein)
# regression_results5 <- merge(regression_results5, proteinmatch, by.x = "Protein", by.y = "idsoma", all.x = T, all.y = F, sort = F)
# 
# ## Save R environment to load in
# save.image(file = "U:/Environments/LMM_environment_soma2.RData")
# load(file = "U:/Environments/LMM_environment_soma2.RData")
# 
# ## Load in environment from Script_linear_mixed models to run code below
# ## Merge LMMs in full list
# compare_lmms <- merge(prot_time_results, prot_time_results_soma2, by = 'Protein')
# compare_lmms$Group <- ifelse(compare_lmms$P.x < 5.8e-5 & compare_lmms$P.y < 5.8e-5, 'Associated in both models',
#                              ifelse(compare_lmms$P.x < 5.8e-5 & compare_lmms$P.y > 5.8e-5, 'Associated in LMM in soma1',
#                                     ifelse(compare_lmms$P.x > 5.8e-5 & compare_lmms$P.y < 5.8e-5, "Associated in LMM in soma2", "Not associated")))
# 
# 
# ## Compare betas from LMM (prot_time_results and prot_time_results_soma2)
# reg <- summary(lm(compare_lmms$`Intervention*v2_Beta.y`~ compare_lmms$`Intervention*v2_Beta.x`))
# rsquared <- paste('R^2 ==', round(reg$adj.r.squared, 2))
# N <- paste('N ==', length(residuals(reg)))
# Fstat <- paste('F ==', reg$fstatistic[1])
# ggplot(compare_lmms, aes( x=`Intervention*v2_Beta.x`, y=`Intervention*v2_Beta.y`)) +
#   geom_point(aes(color = Group )) +
#   geom_smooth(method = "lm", forumla = y ~ x) +
#   xlab("Effect of time on protein using Soma1 data (SDs)") + ylab("Effect of time on protein using soma2 data (SDs)")+
#   annotate("text", x = -1, y = 1, label = rsquared, hjust = 0, vjust =1, parse = TRUE) +
#   annotate("text", x = -1, y = 0.8, label = N, hjust = 0, vjust =1, parse = TRUE) +
#   annotate("text", x = -1, y = 0.6, label = Fstat, hjust = 0, vjust =1, parse = TRUE)
# 
