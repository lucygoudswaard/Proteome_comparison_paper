#### Merging cleaned clinic data and Olink data, adding the freezer variables
library(readxl)
library(tidyr)
library(psych)
library(dplyr)
#source("http://news.mrdwab.com/install_github.R")
#install_github("mrdwab/SOfun")
library(SOfun)
library(lme4)
library(nlme)
#remotes::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze', r$
library(OlinkAnalyze)
library(ggplot2)
library(tidyverse)
library(gtools)
library(qqman)
#install.packages("lubridate")
library(lubridate)
library(stringr)
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats)
library(ggrepel)
library(nlmr)
library(ggpubr)
#library(moosefun)
library(mgcv)


## Set working directory to where the parameter folder is
## Mount project folder on RDSF as well as data warehouse folder

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory
setwd(working_dir)

### Read in metaboprep output and format to give SampleID in first row followed by protein levels
qcdata <- read.table(paste0(prot_release, "2022-05-19/data/metaboprep_release_2021_10_14/filtered_data/LG_bybandsleeve_2021_10_14_Filtered_metabolite_data.txt"), sep = "\t", header=T)
IDs <- as.matrix(rownames(qcdata))
rownames(qcdata) <- NULL
qcdata <- cbind(IDs, qcdata)
names(qcdata)[1] <- "SampleID"

## Read in basic phenotype variables
phenotypes <- read.csv(paste0(clinical, "2022-02-08/data/BBS-2021-09-06_taunton125_analysis_extract_2021-11-15.csv"))
phenotypes$sample_id <- tolower(phenotypes$sample_id)
phenotypes$sample_id <- gsub("-", "", phenotypes$sample_id )

## Merge phenotypes and qcdata dataframes
mydata<- merge(phenotypes, qcdata, by.x = "sample_id", by.y = "SampleID", all.y = TRUE)

## Read in protein names to identify proteins from Olink IDs
olinknames <- read_excel(paste0(prot_release, "2022-05-19/data/raw/20202265_Timpson_NPX_2021-04-27_V2.xlsx"))
OlinkIDs <- olinknames %>% 
  group_by(OlinkID) %>% 
  filter(Index == min(Index)) %>% 
  distinct
OlinkIDs <- OlinkIDs[!duplicated(OlinkIDs$OlinkID), ]
names <- OlinkIDs[, c("OlinkID", "UniProt", "Assay", "MissingFreq")]
names$OlinkID <- tolower(names$OlinkID)
names <- as.data.frame(names)

## LOD info
hist(names$MissingFreq)
mean(names$MissingFreq)
max(names$MissingFreq)
min(names$MissingFreq)
IQR(names$MissingFreq)

## Make study ID factor (this ID identifies individual, combine with time point info to determine
## for each individual which is baseline and 36 months)
mydata <- mydata[order(mydata[, "Study.Identifier"]),]
mydata$Study.Identifier <- as.factor(mydata$Study.Identifier)

## Change visit identifier to numeric
w <- which(mydata$Visit.Identifier == "36 months")
mydata$Visit.Identifier[w] <- "36_months"
#mydata$Visit.Identifier <- as.factor(mydata$Visit.Identifier)
#mydata$Visit.Identifier <- as.numeric(mydata$Visit.Identifier)

## PlateID as factor
str(mydata$olink_Plate.ID)
mydata$olink_Plate.ID <- as.factor(mydata$olink_Plate.ID)

## Set dates for thawing
mydata$baseline_thaw1 <- as.Date("2021-02-08")
mydata$baseline_thaw2 <- as.Date("2021-04-01")

## Change to dates
mydata$Date.research.bloods.taken <- as.Date(mydata$Date.research.bloods.taken)
mydata$Date.bloods.taken <- as.Date(mydata$Date.bloods.taken)
mydata$Randomisation.date <- as.Date(mydata$Randomisation.date)
mydata$Visit.Date <- as.Date(mydata$Visit.Date)

## Create variable for time to first thaw for baseline and endpoint samples
mydata$freezer_time_baseline_thaw1 <- as.numeric(abs(mydata$Date.research.bloods.taken - mydata$baseline_thaw1))
mydata$freezer_time_endpoint_thaw1 <- as.numeric(abs(mydata$Date.bloods.taken - mydata$baseline_thaw1))

## Create variable for time to second thaw
mydata$freezer_time_baseline_thaw2 <- as.numeric(abs(mydata$Date.research.bloods.taken - mydata$baseline_thaw2))
mydata$freezer_time_endpoint_thaw2 <- as.numeric(abs(mydata$Date.bloods.taken - mydata$baseline_thaw2))

## Time from first appointment to second
mydata$days_between_appointments <- as.numeric(abs(mydata$Date.research.bloods.taken - mydata$Date.bloods.taken))
hist(mydata$days_between_appointments)

## Keep IDs mentioned twice
length(unique(mydata[,"Study.Identifier"]))
w <- which(mydata$Study.Identifier == "TAU0244")
mydata <- mydata[-w,]
w <- which(mydata$Study.Identifier == "TAU0234")
mydata <- mydata[-w,]

## Remove IDs that have had no operation
w <- which(mydata$Study.Identifier == "TAU0064" | mydata$Study.Identifier == "TAU0315" | mydata$Study.Identifier == "TAU0920" |  mydata$Study.Identifier == "TAU0072")
mydata <- mydata[-w,]

summary(lm(mydata$oid20325 ~ mydata$Patient.Weight))

### Read in clinical vars
clinic_full <- read.table(paste0(data_intermediate_dir, "Clinical_vars_merged.txt"), sep = "\t", header=T)
#clinic_full$ID_date <- paste0(clinic_full$studyid, "_", clinic_full$visit_date_Baseline)
  
mydata$ID_date <- paste0(mydata$Study.Identifier, "_", mydata$Visit.Date)

### Merge clinic_full and mydata for the wide dataset
all <- merge(clinic_full, mydata, by.x = "studyid", by.y="Study.Identifier", all.y=T) #by = "ID_date"
summary(lm(all$oid20325 ~ all$Patient.Weight))

## Change structure of covariables
all$sex <- as.factor(all$sex)
all$smoke <- as.factor(all$smoke)
all$age_at_rand <- as.numeric(all$age_at_rand)
all$income0 <- as.factor(all$income0)

## Create season variable - summer (1), spring (2), autumn (3), summer (4)
all$Visit_month <- as.numeric(str_replace(all$Visit.Date, pattern = ".*-([^-]*)-.*", replacement = "\\1"))
all$season <- ifelse(all$Visit_month == 12 | all$Visit_month == 1 | all$Visit_month == 2, 1, 
                            ifelse(all$Visit_month == 3 | all$Visit_month == 4 | all$Visit_month == 5, 2,
                                  ifelse(all$Visit_month == 6 | all$Visit_month == 7 | all$Visit_month == 8, 4, 3)))

# Step 1: Inverse normal rank transformation of olink proteins.
#### oid columns are too many
cols <-grepl('oid2', colnames(all))
cols <- which(cols == TRUE)
rnt = as.data.frame(matrix(data = NA, nrow = 238, ncol = 1472 ) )
for (i in 221:1692) {
  rnt[,i-220] <- qnorm((rank(all[,i],na.last="keep")-0.5)/sum(!is.na(all[,i])))
  colnames(rnt)[i-220] <- paste0("rnt_", colnames(all)[i])
  rownames(rnt) <-  all$sample_id
}

# Step 2: fit the ranked data to covariates in a linear model - SHOULD CATEGORICAL VARIABLES BE NUMERIC??
# make dataframe to store residuals from linear regression.
 lm_residuals = as.data.frame(matrix(data = NA, nrow = 238, ncol = 1472 ) )
 # Loop to run linear model for all rank transformed olink proteins and age, sex, plate and season and extract residuals.
  for (i in 1:1472) {
   fit = lm(rnt[,i] ~  all$age_at_rand + all$sex, na.action = na.exclude)
   lm_residuals[,i] <- resid(fit)
  colnames(lm_residuals)[i] <- colnames(rnt)[i]
   rownames(lm_residuals) <-  all$sample_id
 }

# Test normality with Shapiro-Wilk test and extract w-statistics to ensure modelling was appropriate.
shap <-  as.data.frame(matrix(data = NA, nrow = 1472, ncol = 1 ) )
for (i in 1:1472){
  test <- shapiro.test(lm_residuals[,i])
  shap[i,] <- test$statistic
  rownames(shap) <- colnames(lm_residuals)
}
shap <- as.numeric(shap)
hist(shap[,1])

# Combine protein data back with sample info etc
new_data <- merge(all, lm_residuals, by.x="sample_id", by.y="row.names")
summary(lm(new_data$rnt_oid20325 ~ new_data$Patient.Weight))

## Separate protein columns out into timepoints
#small_data <- new_data[ , c(1, 2, 3, 197, 1704:3175)]
new_data$studyid <- str_extract(new_data$ID_date, "[^_]+")
cols <- new_data %>% dplyr::select(grep("rnt_", colnames(new_data)))
colnames <- colnames(cols)
wide_adjusted_olink <- pivot_wider(
  new_data,
  id_cols = studyid,
  names_from = Visit.Identifier,
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = all_of(colnames),
  values_fill = NA,
  values_fn = NULL
)

## Merge these timepoints back to clinic data - wide version
full_data <- merge( clinic_full, wide_adjusted_olink, by = "studyid", all.x=F, all.y = T) 

## Create deltas
cols <-grepl('rnt_', colnames(full_data))
w <- which(cols == TRUE)
start <- w[2]
end <- w[2944]
cols_for_delta <- seq(start, end, 2) ## 2nd entry is baseline therefore i-1 is 36 months

olink_deltas <- as.data.frame(matrix(nrow=119, ncol=1472))
for (i in cols_for_delta) {
  colnames(olink_deltas)[which(cols_for_delta == i)] <- colnames(full_data)[i]
  olink_deltas[,which(cols_for_delta == i)] <- full_data[,i-1] - full_data[,i]
  colnames(olink_deltas) <- gsub("rnt_", "", colnames(olink_deltas))
  colnames(olink_deltas) <- gsub("_Baseline", "_change", colnames(olink_deltas))
  rownames(olink_deltas) <- full_data$studyid
}

data_large <- merge(full_data, olink_deltas, by.x = "studyid", by.y = "row.names")

## Mock regression - using IGFBP2 (oid20325)
lm <- lm(data_large$oid20325_change ~ data_large$BMI_change + data_large$age_at_rand + data_large$sex)
summary(lm)
t.test(data_large$rnt_oid20325_Baseline, data_large$rnt_oid20325_36_months, paired = T )

## Plot BMI change and IGFBP2 change (oid20325)
ggplot(data_large, aes(x = BMI_change, y = oid20325_change)) +
  geom_point(aes(x = BMI_change, y = oid20325_change)) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  xlab("BMI change (kg/m2)") + ylab("IGFBP2 change (SDs)") + ggtitle("") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) 

## Same with leptin (oid20187)
lm <- lm(data_large$oid20187_change ~ data_large$BMI_change + data_large$age_at_rand + data_large$sex)
summary(lm)

ggplot(data_large, aes(x = BMI_change, y = oid20187_change)) +
  geom_point(aes(x = BMI_change, y = oid20187_change)) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  xlab("BMI change (kg/m2)") + ylab("Leptin change (SDs)") + ggtitle("") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) 

lm <- lm(data_large$rnt_oid20187_36_months ~ data_large$bmi_kg_m2_36_months + data_large$rnt_oid20187_Baseline)
summary(lm)

### Run regression over all protein change cols
linear_model <- function(wdata, dependent, independent, covariables) {
  ##############################
  ### 1. Define Model Data Frame
  ##############################
  if(is.na(covariables[1])){
    model_variables = c( dependent, independent )
  } else {
    model_variables = c( dependent, covariables, independent)  
  }
  
  mod_data = wdata[, c(model_variables)]
  ##############################
  ### 2. Perform Linear Model 
  ##############################
  if( is.na( covariables[1] ) ){
    form = formula(paste0(dependent ," ~ ", independent ))  
  } else {
    form = formula(paste0(dependent ," ~ ", independent ," + ", paste0(covariables, collapse = " + ")   ))  
  }
  
  lm_mod <- lm(form, data = mod_data)
  
  #################
  ## 3. summary stats
  #################
  s = summary(lm_mod)
  ## sample size
  n = length(residuals(s)); names(n) = "n_lm"
  ## Adjusted R-squared, to compare between models with a different number of predictors
  rsq = s$adj.r.squared; names(rsq) = "rsq_adj_lm"
  ## Dependent Variable effect estimates
  beta = s$coefficients[2,1]; names(beta) = "beta_lm"
  se = s$coefficients[2,2]; names(se) = "se_lm"
  pval = s$coefficients[2,4]; names(pval) = "P_lm"
  lm_results = c(n, rsq, beta, se, pval)
  return(lm_results)
}

## Run all protein columns baseline
cols <-grepl('_change', colnames(data_large))
w <- which(cols == TRUE)
traits = colnames(data_large)[w]
traits <- traits[6:1477]

## Regression for BMI change and protein change
#data_large$sex <- as.numeric(data_large$sex)
bmi_change_protein_change = t( sapply(traits, function(trait){
  linear_model(wdata = data_large, 
               dependent = trait, 
               independent = "BMI_change",
               covariables =  c("sex", "age_at_rand"))
}) )

row.names(bmi_change_protein_change) <- sub("_.*", "", row.names(bmi_change_protein_change) )
bmi_change_protein_change_full <- merge(names, bmi_change_protein_change,  by.x='OlinkID', by.y="row.names")
bmi_change_protein_change_full <- bmi_change_protein_change_full[order(bmi_change_protein_change_full[, "P_lm"]),]

### Write results
#write.table(bmi_change_protein_change_full, file = paste0(data_output_dir, "tables/06_BMI_change_prot_change_reg.txt"), sep = "\t", col.names=T, row.names=F)

### Volcano
bmi_change_protein_change_full$Category <- ifelse(bmi_change_protein_change_full$P_lm <6.2e-5, "Associated", "Not associated")
volcano <- ggplot(data=bmi_change_protein_change_full, 
                  aes(x=beta_lm, y=-1*log10(P_lm), col=Category, 
                      label=ifelse(P_lm < 6.2e-5, Assay, ""))) + 
  geom_text(size =2) +
  geom_point(size =0.5) +
  theme_minimal() + 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=-log10(6.2e-5), col="red") + 
  xlab("Beta coefficient (SD change per kg/m2 increase in BMI)") +
  ylab(expression("-log"[10]*"(p-value)"))

volcano

#### Linear mixed model
### Data in long format - use new_data
LMM_data <- new_data

## Make Visit identifier as factor. 1 = 36 months, Baseline = 2
w <- which(LMM_data$Visit.Identifier == "Baseline")
LMM_data$Visit_numeric[w] <- 1 
w <- which(LMM_data$Visit.Identifier == "36_months")
LMM_data$Visit_numeric[w] <- 2 

## Run linear mixed model - example with leptin (oid20187)
### LMER
LMM <- lme4::lmer(LMM_data$rnt_oid20187 ~ LMM_data$Visit_numeric + (1|LMM_data$studyid))
summary(LMM)

## IGFBP2 oid20325
LMM <- lme4::lmer(LMM_data$rnt_oid20325 ~ LMM_data$Visit_numeric + (1|LMM_data$studyid))
sum <- summary(LMM)

## Run for all proteins
## Need to remove any participants who have data only at one timepoint
cols <-grepl('rnt_', colnames(LMM_data))
w <- which(cols == TRUE)
x <- as.data.frame(matrix(ncol=2, nrow=1472))
for (i in w) {
  x[which(w==i),1] <- colnames(LMM_data)[i]
  x[which(w==i),2] <- sum(is.na(LMM_data[,i]))
}
table(x[,2])

w <- which(is.na(LMM_data$rnt_oid21502))
LMM_data$studyid[w]

## rows 47 and 194 are TAU0822, turn to NA for this protein
#LMM_data$rnt_oid21502[47] <- NA 
#LMM_data$rnt_oid21502[190] <- NA
LMM_data <- LMM_data[-c(47, 188),]

lmer_results <- as.data.frame(matrix(ncol=6, nrow=1472))
colnames(lmer_results) <- c("Protein", "Fixed_effects_estimate", "SE", "t", "N", "P")
cols <-grepl('rnt_oid', colnames(LMM_data))
w <- which(cols == TRUE)

for (i in w){
  pt <- lme4::lmer(LMM_data[,i] ~ LMM_data$Visit_numeric + (1|LMM_data$studyid))
  protein_time <- summary(pt)
  lmer_results[,1] <- colnames(LMM_data)[w]
  lmer_results[i-(w[1]-1),2] <- protein_time$coefficients[2,1]
  lmer_results[i-(w[1]-1),3] <- protein_time$coefficients[2,2]
  lmer_results[i-(w[1]-1),4] <- protein_time$coefficients[2,3]
  lmer_results[i-(w[1]-1),5] <- length(residuals(pt))
  fitA <- lmer(LMM_data[,i] ~ LMM_data$Visit_numeric + (1|LMM_data$studyid), REML=FALSE)
  fitB <- lmer(LMM_data[,i] ~ (1|LMM_data$studyid), REML=FALSE)
  compare <- anova(fitA,fitB) 
  ### Chisq is [2,6], pvalue is [2,8] of compare
  lmer_results[i-(w[1]-1),6] <- compare[2,8]
}

lmer_results$Protein2 <- sub(".*_", "", lmer_results$Protein )
lmer_results_full <- merge(names, lmer_results,  by.x='OlinkID', by.y="Protein2")
lmer_results_full <- lmer_results_full[order(lmer_results_full[, "P"]),]

## LOD of associated results
w <- which(lmer_results_full$P < 6.2e-5)
assoc <- lmer_results_full[w,]
mean(assoc$MissingFreq)
min(assoc$MissingFreq)
max(assoc$MissingFreq)

#write.table(lmer_results_full, file = paste0(data_output_dir, "tables/06_time_rnt_protein_LMM.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")


## Volcano plot for linear model
lmer_results_full$Category <- ifelse(lmer_results_full$P <6.2e-5, "Associated", "Not associated")
volcano <- ggplot(data=lmer_results_full, 
                  aes(x=Fixed_effects_estimate, y=-1*log10(P), col=Category)) + 
                     # label=ifelse(P < 1e-10, Assay, ""))) + 
  geom_vline(xintercept=0, col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(6.2e-5), col="grey", linetype= "dashed") + 
  geom_label_repel(aes(label = ifelse(P < 6.2e-5, Assay, "")), max.overlaps = 20, size = 1.5) +
 # geom_text(size =2) +
  geom_point(size =0.5) +
  theme_minimal() + 

  xlab("By-Band-Sleeve intervention effect on protein levels (SDs)") +
  ylab(expression("-log"[10]*"(p-value)"))


#pdf(file = paste0(data_output_dir, "figures/06_volcano_BBS_LMM_boxed_labels.pdf"), width = 6, height = 6)
volcano
dev.off()

### Try linear mixed model adjusting for age and sex - how does this compare without these in the model?
lmer_results2 <- as.data.frame(matrix(ncol=6, nrow=1472))
colnames(lmer_results2) <- c("Protein", "Fixed_effects_estimate", "SE", "t", "N", "P")
cols <-grepl('rnt_oid', colnames(LMM_data))
w <- which(cols == TRUE)

for (i in w){
  pt2 <- lme4::lmer(LMM_data[,i] ~ LMM_data$Visit_numeric + (1|LMM_data$studyid) + LMM_data$age_at_rand + LMM_data$sex)
  protein_time2 <- summary(pt2)
  lmer_results2[,1] <- colnames(LMM_data)[w]
  lmer_results2[i-(w[1]-1),2] <- protein_time2$coefficients[2,1]
  lmer_results2[i-(w[1]-1),3] <- protein_time2$coefficients[2,2]
  lmer_results2[i-(w[1]-1),4] <- protein_time2$coefficients[2,3]
  lmer_results2[i-(w[1]-1),5] <- length(residuals(pt2))
  fitA <- lmer(LMM_data[,i] ~ LMM_data$Visit_numeric + (1|LMM_data$studyid) + LMM_data$age_at_rand + LMM_data$sex, REML=FALSE)
  fitB <- lmer(LMM_data[,i] ~ (1|LMM_data$studyid) + LMM_data$age_at_rand + LMM_data$sex, REML=FALSE)
  compare2 <- anova(fitA,fitB) 
  ### Chisq is [2,6], pvalue is [2,8] of compare
  lmer_results2[i-(w[1]-1),6] <- compare2[2,8]
}

lmer_results2$Protein2 <- sub(".*_", "", lmer_results2$Protein )
lmer_results2_full <- merge(names, lmer_results2,  by.x='OlinkID', by.y="Protein2")
lmer_results2_full <- lmer_results2_full[order(lmer_results2_full[, "P"]),]

#write.table(lmer_results2_full, file = paste0(data_output_dir, "tables/06_time_rnt_protein_LMM_age_sex.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")

## Merge results
lmer_combined <- merge(lmer_results_full, lmer_results2_full, by = "OlinkID")

### Plot estimates and SEs against each other
lm<- summary(lm(data=lmer_combined, Fixed_effects_estimate.y ~ Fixed_effects_estimate.x))
#pdf(file = "figures/06_lmer_method_comparison.pdf", width = 12, height = 10)
ggplot(lmer_combined, aes( x = Fixed_effects_estimate.x, y = Fixed_effects_estimate.y )) +
  geom_point(aes( size = 1)) +
  geom_errorbar(aes(ymin = Fixed_effects_estimate.y-SE.y,ymax = Fixed_effects_estimate.y+SE.y), alpha = 0.1) + 
  geom_errorbarh(aes(xmin = Fixed_effects_estimate.x-SE.x,xmax = Fixed_effects_estimate.x+SE.x), alpha = 0.1) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  xlab("BBS linear mixed model estimate (SD)") + ylab("BBS linear mixed model estimate with age and sex (SD)") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) +
  annotate("text", x = -0.8, y = 1, label = "N=1472", size = 5) +
  annotate("text", x = -0.8, y = 0.9, label = ("R squared = 1"), size = 5) +
  annotate("text", x = -0.8, y = 0.8, label = ("p<2.2e-16"), size = 5)
#dev.off()
