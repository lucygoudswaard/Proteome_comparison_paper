######## Explore the clinic variables in BBS #########
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
library(gtsummary)
library(stringr)

## Set working directory to where the parameter folder is

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory
setwd(working_dir)

### Read in clinic files
baseline <- read.csv(paste0(clinical, "2022-02-08/data/20220118_database_extract/metabolomics_TAUsample_baseline_18JAN2022.csv"), sep = "|", na.strings = ".")
bloods <- read.csv(paste0(clinical, "2022-02-08/data/20220118_database_extract/metabolomics_TAUsample_bloods_18JAN2022.csv"), sep = "|", na.strings = ".")
followup <- read.csv(paste0(clinical, "2022-02-08/data/20220118_database_extract/metabolomics_TAUsample_followup_18JAN2022.csv"), sep = "|", na.strings = ".")
meds <- read.csv(paste0(clinical, "2022-02-08/data/20220118_database_extract/metabolomics_TAUsample_meds_18JAN2022.csv"), sep = "|", na.strings = ".")
barcodes <- read.csv(paste0(clinical, "2022-02-08/data/20220118_database_extract/metabolomics_TAUsample_barcodes_18JAN2022.csv"), sep = "|", na.strings = ".")

### Remove any duplicated rows
## Drop ELF and ELF date to be able to do this
bloods <- subset(bloods, select = -c(ELF,ELF_date) )
bloods <- unique(bloods)

## Merge unique IDs into bloods
phenotypes <- read.csv(paste0(clinical, "2022-02-08/data/BBS-2021-09-06_taunton125_analysis_extract_2021-11-15.csv"))
phenotypes$sample_id <- tolower(phenotypes$sample_id)
phenotypes$sample_id <- gsub("-", "", phenotypes$sample_id )
w <- which(phenotypes$Visit.Identifier == "36 months")
phenotypes$Visit.Identifier[w] <- "36_months"
phenotypes$ID <- paste0(phenotypes$Study.Identifier, "_", phenotypes$Visit.Identifier)

## Restrict bloods to timepoints of interest
bloods$Timepoint <- "Baseline"
w <- which(bloods$VisitID == 6)
bloods$Timepoint[w] <- "36_months"
bloods$ID <- paste0(bloods$StudyID, "_", bloods$Timepoint)

### merge timepoints into blood
bloods <- merge(phenotypes, bloods, by="ID")

## Change blood into wide format
bloodswide <- pivot_wider(
  bloods,
  id_cols = c(StudyID),
  names_from = Timepoint,
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c( hb, plate, wbc, rbc, hct, mch, mcv, mchc, lymph, neutro, sodium, potass,
                  urea, creat, alt, alp, albumin, bili, bili_sign, hba1c, serumiron,
                  ferrit, vitb12, folate, folate_sign, hydroxy, calc, parathyroid, crp,
                  mag, phos, totalprot, blood_date, blood_time, chol, hdlc, ldlc, triglyc,
                  fastgluc, fblood_date, fblood_time, blood_liver, blood_liver_date, blood_liver_time,
                  ELF_Fibrosis, blood_research, blood_research_date, blood_research_time),
  values_fill = NA,
  values_fn = NULL
)

## Restrict follow up to timepoint 0 and 6
w <- which(followup$visitid == 0 | followup$visitid == 6)
followup <- followup[w,]

followup$Timepoint <- "Baseline"
w <- which(followup$visitid == 6)
followup$Timepoint[w] <- "36_months"

followup_wide <- pivot_wider(
  followup,
  id_cols = studyid,
  names_from = Timepoint,
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c(visit_date, wgt, preop_wgt_date, waist_1, waist_2, waist_3,
                  waist_4, sbp_1, dbp_1, sbp_2, dbp_3, sbp_3, dbp_3),
  values_fill = NA,
  values_fn = NULL
)

### Restrict meds timepoints and make meds in a wider format
meds$Timepoint <- "Baseline"
w <- which(meds$visitid == 6)
meds$Timepoint[w] <- "36_months"

meds_wide <- pivot_wider(
  meds,
  id_cols = StudyID,
  names_from = c(Timepoint, drug_cat),
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c(drug_name),
  #values_from = c(drug_name_comp, visitdate, drug_id, drug_cat, drug_dose, drug_freq,
                  #drug_freqoth, drug_name, drug_record_type, drug_date, recoded, drug_cat_old, 
                  #sub_type1),
  values_fill = NA,
  values_fn = NULL
)

cols <- 2:23
for (i in cols) {
  meds_wide[,i] <- ifelse(meds_wide[,i]=="NULL", "No", "Yes")
}

## Merge all data
mydata_prep <- merge(baseline, followup_wide, by="studyid")
mydata_prep2 <- merge(mydata_prep, bloodswide, by.x="studyid", by.y="StudyID")
mydata <- merge(mydata_prep2, meds_wide, by.x="studyid", by.y ="StudyID")

## Derive baseline and endpoint BMI
mydata$hgt_m <- mydata$hgt/100
mydata$bmi_kg_m2_Baseline <- mydata$wgt_Baseline / ((mydata$hgt_m)*(mydata$hgt_m)) 
mydata$bmi_kg_m2_36_months <- mydata$wgt_36_months / ((mydata$hgt_m )*(mydata$hgt_m)) 

## BMI change
mydata$BMI_change <- mydata$bmi_kg_m2_36_months - mydata$bmi_kg_m2_Baseline

## Wight change
mydata$weight_change <- mydata$wgt_36_months - mydata$wgt_Baseline

## Waist circumference averages
mydata$waist_3_Baseline <- as.numeric(mydata$waist_3_Baseline)
mydata$waist_4_Baseline <- as.numeric(mydata$waist_4_Baseline)
mydata$waist_average_Baseline <- rowMeans(mydata[, c("waist_1_Baseline", 
    "waist_2_Baseline", "waist_3_Baseline", "waist_4_Baseline")], na.rm=TRUE ) ## No missing

mydata$waist_3_36_months <- as.numeric(mydata$waist_3_36_months)
mydata$waist_4_36_months <- as.numeric(mydata$waist_4_36_months)
mydata$waist_average_36_months <- rowMeans(mydata[, c("waist_1_36_months", 
    "waist_2_36_months", "waist_3_36_months", "waist_4_36_months")], na.rm=TRUE ) ## 1 missing

## Derive waist circumference change??
mydata$waist_change <- mydata$waist_average_36_months - mydata$waist_average_Baseline

## Blood pressure averages. Change BP variables to continuous.
voi <- mydata[,c("sbp_1_Baseline", "sbp_1_36_months", "dbp_1_Baseline", "dbp_1_36_months",
               "sbp_2_Baseline", "sbp_2_36_months", "dbp_3_Baseline", "dbp_3_36_months",
               "sbp_3_Baseline", "sbp_3_36_months")]

for (i in 1:10) {
  voi[,i] <- as.numeric(voi[,i])
}

## Derive average blood pressures from repeats
mydata$sbp_average_baseline <- rowMeans(voi[, c("sbp_1_Baseline", 
    "sbp_2_Baseline", "sbp_2_Baseline")], na.rm=TRUE )

mydata$sbp_average_36_months <- rowMeans(voi[, c("sbp_1_36_months", 
    "sbp_2_36_months", "sbp_2_36_months")], na.rm=TRUE ) ## 1 NA

mydata$dbp_average_baseline <- rowMeans(voi[, c("dbp_1_Baseline", 
  "dbp_3_Baseline")], na.rm=TRUE )

mydata$dbp_average_36_months <- rowMeans(voi[, c("dbp_1_36_months", 
  "dbp_3_36_months")], na.rm=TRUE ) ## 1 NA

## Derive BP change
mydata$sbp_change <- mydata$sbp_average_36_months - mydata$sbp_average_baseline 
mydata$dbp_change <- mydata$dbp_average_36_months - mydata$dbp_average_baseline

## Derive time with diabetes before enrollment - combine month and year variables. 
## For ease, the start of each month is used as there is no day variable
mydata$diabetes_month_year <- na.omit(paste0("1/", mydata$diab_date_m, "/",  mydata$diab_date_y))
w <- which(mydata$diabetes_month_year == "1/NA/NA")
mydata$diabetes_month_year[w] <- NA
mydata$diabetes_month_year <- as.Date(mydata$diabetes_month_year, format="%d/%m/%Y")
mydata$visit_date_Baseline <- as.Date(mydata$visit_date_Baseline, format="%d/%m/%Y")

mydata$diabetes_duration <- interval(mydata$diabetes_month_year, mydata$visit_date_Baseline) %/% months(1)

## More time variables
# Time from baseline sampling to surgery - baseline visit date to operation date
mydata$opdate <- as.Date(mydata$opdate, format = "%d/%m/%Y")
mydata$baseline_to_op <- difftime(mydata$opdate, mydata$visit_date_Baseline)
mydata$baseline_to_op <- as.numeric(mydata$baseline_to_op)
min(mydata$baseline_to_op, na.rm = T) #18
max(mydata$baseline_to_op, na.rm = T) #1304!

# Time from operation to end sampling - op date to visit 6
mydata$visit_date_36_months <- as.Date(mydata$visit_date_36_months, format="%d/%m/%Y")
mydata$op_to_end <- difftime(mydata$visit_date_36_months, mydata$opdate)
mydata$op_to_end <- as.numeric(mydata$op_to_end)
min(mydata$op_to_end, na.rm = T) #-252
max(mydata$op_to_end, na.rm = T) #1335

# Time from baseline sampling to end sampling - visit 0 to visit 6
mydata$baseline_to_end <- difftime(mydata$visit_date_36_months, mydata$visit_date_Baseline)
mydata$baseline_to_end <- as.numeric(mydata$baseline_to_end)
min(mydata$baseline_to_end, na.rm = T) #1043
max(mydata$baseline_to_end, na.rm = T) #1426

# Sample time in freezer (derived in other script)
mydata$date_thawed <- as.Date("2021-02-08")
mydata$blood_research_date_Baseline <- as.Date(mydata$blood_research_date_Baseline, format = "%d/%m/%Y")
mydata$time_in_freezer_baseline <- difftime(mydata$date_thawed, mydata$blood_research_date_Baseline)
mydata$fblood_date_36_months <- as.Date(mydata$fblood_date_36_months, format = "%d/%m/%Y")
mydata$time_in_freezer_endpoint <- difftime(mydata$date_thawed, mydata$fblood_date_36_months)
mydata$time_in_freezer_baseline <- as.numeric(mydata$time_in_freezer_baseline)
mydata$time_in_freezer_endpoint <- as.numeric(mydata$time_in_freezer_endpoint)
min(mydata$time_in_freezer_baseline, na.rm = T) #1200
max(mydata$time_in_freezer_baseline, na.rm = T) #2905
min(mydata$time_in_freezer_endpoint, na.rm = T) #95
max(mydata$time_in_freezer_endpoint, na.rm = T) #1803

# Time difference between weight measures and sampling - visit date to blood research date (should mostly be the same)
mydata$visit_date_Baseline %in% mydata$blood_research_date_Baseline ## FALSE for two ppts
mydata$visit_date_36_months %in% mydata$fblood_date_36_months ## 12 FALSE - probably for COVID reasons
mydata$visit_date_to_bloods_baseline <- difftime(mydata$blood_research_date_Baseline, mydata$visit_date_Baseline, units = "days")
mydata$visit_date_to_bloods_endpoint <- difftime(mydata$fblood_date_36_months, mydata$visit_date_36_months, units = "days")

### Write clinical vars as intermediate dataset
#write.table(mydata, file = paste0(data_intermediate_dir, "Clinical_vars_merged.txt"), sep = "\t", col.names = T, row.names = F )

## Subset all continuous variables and remove those that are just NAs
data_hist <- as.data.frame(lapply(mydata,as.numeric))
data_hist <- data_hist[,colSums(is.na(data_hist))<nrow(data_hist)]

## Remove all the categorical variables from histograms
data_hist <- data_hist[, -which(names(data_hist) %in% c ("sex", "ethnicity", "employ", "income0", "income_36m", "type1_diab0", "type2_diab0",
  "diab_date_m", "angina0", "prevmi0", 'prevcabg0', "prevstroke0", "prevclaud0", "nyha0",
  "smoke", "visit_date_Baseline", "folate_sign_36_months", "diabetes_month_year"))]

#pdf(file = paste0(data_output_dir, "figures/05_Clinic_vars_histograms.pdf"), height = 10, width = 7)
par(mfrow=c(3,2))
for (i in 1:98) {
  DataDescribed = psych::describe(data_hist[,i]) 
  meanvar<-DataDescribed$mean
  medianvar<-DataDescribed$median
  minvar<-DataDescribed$min
  maxvar<-DataDescribed$max
  kurtosisvar<-DataDescribed$kurtosis
  skewnessvar<-DataDescribed$skew
  N<- nrow(data_hist) - sum(is.na(data_hist[,i]))
  missingness <- (sum(is.na(data_hist[,i]))/nrow(data_hist))*100
  
  a<-density(data_hist[,i], na.rm=T)
  thresholdx<-(maxvar+(maxvar/100))
  thresholdy<-min(a$y)+(max(a$y)/4)
  
  hist(data_hist[,i], col="red",main=(names(data_hist)[i]),prob=TRUE,xlab="peak area") 
  lines(density(data_hist[,i], na.rm = TRUE),col="blue", lwd=2)
  text(thresholdx,thresholdy, cex=0.6, 
       paste("N=", N, "\npercent missing=", 
             signif(missingness, 3), "\nmin=", 
             signif(minvar, 3), " \nmax=",
             signif(maxvar, 3), 
             "\nmean=", 
             signif(meanvar, 3), " \nmedian=", signif(medianvar, 3), 
             "\nkurt=", 
             signif(kurtosisvar, 3), " \nskew=", 
             signif(skewnessvar, 3), sep = ''), pos = 3,xpd = NA)
}
dev.off() 

## Boxplot of BMI
bmi_box_data <- as.data.frame(matrix(nrow=250 , ncol=2))
colnames(bmi_box_data) <- c("BMI", "Timepoint")
bmi_box_data[,1] <- c(mydata$bmi_kg_m2_Baseline, mydata$bmi_kg_m2_36_months)
bmi_box_data[1:125,2] <- "Baseline"
bmi_box_data[126:250, 2] <-  "Endpoint (36m)"

colors <- c("red", "pink")
#pdf(file = paste0(data_output_dir, "figures/05_BMI_timepoint_histograms.pdf"), height = 4, width = 4)
bmibox <- boxplot(BMI ~ Timepoint, data=bmi_box_data, boxfill = colors,
  xlab = "Timepoint", ylab = expression("BMI (kg/m"^"2"*")")) 
 dev.off()                 
 
## Summary of characteristics of dataset used for main analysis. Remove IDs that are removed from the main analysis.
## IDs have been replaced with numbers so they are not identifiable
cat_characteristics <- mydata
cat_characteristics <- as.data.frame(cat_characteristics)
w <- which(cat_characteristics$studyid == "1" | cat_characteristics$studyid == "2" |
             cat_characteristics$studyid == "3" |  cat_characteristics$studyid == "4" | 
             cat_characteristics$studyid == "5" | cat_characteristics$studyid == "6" | 
             cat_characteristics$studyid == "7")
cat_characteristics <- cat_characteristics[-w,]

## For included participants, what is the mean and rage of operation to endpoint
mean(cat_characteristics$op_to_end / 365.25 * 12, na.rm = T)
min(cat_characteristics$op_to_end / 365.25 * 12, na.rm = T)
max(cat_characteristics$op_to_end / 365.25 * 12, na.rm = T)

cat_characteristics <- cat_characteristics[, c("age_at_rand",  "sex", "bmi_kg_m2_Baseline", "sbp_average_baseline", "dbp_average_baseline", "smoke", "diabetes", "diabetes_duration",
                                               "ethnicity", "employ", "income0", "type1_diab0", "type2_diab0", "triglyc_Baseline", "ldlc_Baseline", "hdlc_Baseline", "chol_Baseline" )]
cat_characteristics$diabetes <- as.factor(cat_characteristics$diabetes)
cat_characteristics$triglyc_Baseline <- as.numeric(cat_characteristics$triglyc_Baseline)
cat_characteristics$ldlc_Baseline <- as.numeric(cat_characteristics$ldlc_Baseline)
cat_characteristics$hdlc_Baseline <- as.numeric(cat_characteristics$hdlc_Baseline)
cat_characteristics$chol_Baseline <- as.numeric(cat_characteristics$chol_Baseline)  

### Categorical table
table1 <- tbl_summary(cat_characteristics, statistic = list(all_continuous() ~ "{mean} ({sd})",
                                                                all_categorical() ~ "{n} / {N} ({p}%)"), missing = "no")

table1 %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(data_output_dir, "figures/05_characteristic_baseline.html"))

table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = paste0(data_output_dir, "figures/05_characteristic_baseline.docx"))

### Categorical table by sex
table2 <- tbl_summary(cat_characteristics, statistic = list(all_continuous() ~ "{mean} ({sd})",
                                                            all_categorical() ~ "{n} / {N} ({p}%)"), missing = "no", by="sex") 
table2 <- table2 %>%   add_p(pvalue_fun = function(x) signif(x, 2))
 
## write to file
table2 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = paste0(data_output_dir, "figures/05_characteristic_baseline_by_sex.docx"))

## Same table but with categories collapsed
data_wide_collapsed <- cat_characteristics
data_wide_collapsed$ethnicity <- as.factor(data_wide_collapsed$ethnicity)
data_wide_collapsed$ethnicity <- fct_collapse(data_wide_collapsed$ethnicity,
                                              "1" = "1",
                                              "2"= c("2", "3", "4", "5"))
data_wide_collapsed$employ <- as.factor(data_wide_collapsed$employ)
data_wide_collapsed$employ <- fct_collapse(data_wide_collapsed$employ,
                                                        "1" = "1",
                                                        "2" = "2",
                                                        "3" = "3",
                                                        Unemployed = c('4', '5', '6', '7', '8'))
data_wide_collapsed$income0 <- as.factor(data_wide_collapsed$income0)
data_wide_collapsed$income0 <- fct_collapse(data_wide_collapsed$income0,
                                           'Low' = c("1", "2"),
                                           'Medium' = c("3", "4"),
                                           "High" = c( '5', '6', '7', '8'))

data_wide_collapsed$smoke <- as.factor(data_wide_collapsed$smoke)
data_wide_collapsed$smoke <- fct_collapse(data_wide_collapsed$smoke,
                                                     '2' = c("2", "3"))

data_wide_collapsed$type1_diab0 <- as.factor(data_wide_collapsed$type1_diab0)
data_wide_collapsed$type2_diab0 <- as.factor(data_wide_collapsed$type2_diab0)

table3 <- tbl_summary(data_wide_collapsed, by = sex, missing = "no",
                      statistic = list(all_continuous() ~ "{mean} ({sd})",
                                       all_categorical() ~ "{n} ({p}%)")) %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) signif(x, 2))

# write to file
table3 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = paste0(data_output_dir, "figures/05_characteristic_baseline_by_sex_collapsed_cats.docx"))

