####### Comparing three studies - INTERVAL, DiRECT, BBS, using LMM in all DiRECT data  5/10/22

########### DiRECT and INTERVAL estimate omparison
library(zip)
#install.packages("openxlsx")
library(openxlsx)
#install.packages("data.table")
library(data.table)
#install.packages("systemfit")
library(systemfit)
#devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
#install.packages("devtools")
library(devtools)
#install.packages("dplyr")
library(dplyr)
#install.packages("cowplot")
library(cowplot)
#install.packages("reshape2")
library(reshape2)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
#install.packages("tidyr")
library(tidyr)
#install.packages("purrr")
library(purrr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("broom")
library(broom)
#install.packages("generics")
library(generics)
#devtools::install_github("tidyverse/glue")
library(glue)
#devtools::install_github(repo = "amices/mice")
library(utils)
library(mice)
#install.packages("purrr")
library(purrr)
#install.packages("zoo")
library(zoo)
#install.packages("psych")
library(psych)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("doppelgangR")
library(doppelgangR)
# install.packages("gtsummary")
library(gtsummary)
library(corrplot)
library(ggforestplot)
library(dplyr)
#install.packages("deming")
library(deming)
library(qqman)
library(rgl)
library(scatterplot3d)
library(readxl)
library(ggrepel)

setwd("file/path")

DiRECT <- read_excel("bmi_change_protein_main_delta_rnt_26.4.21.xlsx", sheet = "2SLS_BMI_protein")
DiRECT <- DiRECT[-c(16)]
DiRECTobs <- read_excel("bmi_change_protein_main_delta_rnt_26.4.21.xlsx", sheet = "Obs_BMI_change_protein")
INTERVALobs <- read_excel("INETRVAL_MR_obs_for_comparison.xlsx", sheet = "Obs")
INTERVALMR <- read_excel("INETRVAL_MR_obs_for_comparison.xlsx", sheet = "MR")

INTERVAL <- merge(INTERVALobs, INTERVALMR, by = "ID")
INTERVAL <- distinct(INTERVAL, Target.x, .keep_all = T)

## Restrict to unique Target names
DiRECT <- distinct(DiRECT,Target, .keep_all = T) 
DandI <- merge(DiRECT, INTERVAL, by.x = "Target", by.y = "Target.x", all= FALSE)

### Read in some Olink results from ByBandSleeve to merge with the Soma data
olink <- read.table(file = "file/path/06_BMI_change_prot_change_reg.txt", sep = "\t", header = T)
all_models <- merge(olink, DandI, by.x = "UniProt", by.y = "UniProtID.x", all = F)
all_models$Category <- ifelse(all_models$Pval < 0.001 & all_models$P_val.y < 0.001 & all_models$P_lm < 0.001 , "Associated in all", "Not associated in all")

### Read in alternative LMM files
BBS_LMM <- read.table(file = "file/path/06_time_rnt_protein_LMM.txt", sep = "\t", header = T)
### How many proteins increase with 
w <- which(BBS_LMM$P < 6.2e-5) # 191
Associated_BBS <- BBS_LMM[w,]
Direction <- sign(Associated_BBS$Fixed_effects_estimate)
pos <- length(which(Direction == 1)) ## 118
neg <- length(which(Direction == -1)) ## 73

BBS_LMM <- distinct(BBS_LMM,UniProt, .keep_all = T) 

## Read in linear mixed model for DiRECT
DiRECT_LMM <- read_excel("lmer_protein_time_soma2_30.11.22.xlsx")
w <- which(DiRECT_LMM$P < 2.1e-5) ## 216
Associated_DiRECT <- DiRECT_LMM[w,]
Direction <- sign(Associated_DiRECT$`Intervention*v2_Beta`)
pos <- length(which(Direction == 1)) ## 120
neg <- length(which(Direction == -1)) ## 96

## DiRECT volcano plot
DiRECT_LMM$Category <- ifelse(DiRECT_LMM$P <2.1e-5, "Associated", "Not associated")

volcano <- ggplot(data=DiRECT_LMM, 
                  aes(x=`Intervention*v2_Beta`, y=-1*log10(P), col=Category)) +
  #  label=ifelse(P < 1e-15, Target, ""))) + 
  geom_label_repel(aes(label = ifelse(P < 1e-10, Target, "")), max.overlaps = 20, size = 1.5) +
  geom_vline(xintercept=0, col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(2.1e-5), col="grey", linetype="dashed") + 
  # geom_text(size =2) +
  geom_point(size =0.5) +
  theme_minimal() + 
  xlab("DiRECT intervention effect on protein levels (SDs)") +
  ylab(expression("-log"[10]*"(p-value)"))

#pdf(file = "volcano_DiRECT_box_labs_soma2.pdf", height = 6, width = 6)
volcano
dev.off()

## Add in protein annotation to file
#Annotation <- DiRECT[,c("Protein", "SomaId", "TargetFullName", "Target", "UniProt")]
#DiRECT_LMM <- merge(DiRECT_LMM, Annotation, by = "Protein")
DiRECT_LMM <- distinct(DiRECT_LMM,Target, .keep_all = T) 

## Merge
INTERVAL_DiRECT <- merge(INTERVALMR, DiRECT_LMM, by.x = "UniProtID", by.y = "UniProt", all = FALSE)
INTERVAL_DiRECT <- distinct(INTERVAL_DiRECT,Target.x, .keep_all = T) 

## Merge BBS
all_LMMs <- merge(INTERVAL_DiRECT, BBS_LMM, by.x = "UniProtID", by.y = "UniProt", all = FALSE)

#Comparison <- similar_estimates(weight_loss$`Intervention*v2_Beta`, weight_loss$Fixed_effects_estimate)
weight_loss <- merge(DiRECT_LMM, BBS_LMM, by="UniProt")
weight_loss$Association <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5,  "Associated in both studies", "Not associated in both studies")

weight_loss$Opposite_direction <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` < 0 & weight_loss$Fixed_effects_estimate > 0, "Opposite", "Not opposite")
weight_loss$Group <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` < 0 & weight_loss$Fixed_effects_estimate < 0, "Associated in both studies",
                            ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` > 0 & weight_loss$Fixed_effects_estimate > 0, "Associated in both studies", "Not associated in both studies"))
weight_loss$multiply <- weight_loss$`Intervention*v2_Beta` * weight_loss$Fixed_effects_estimate

weight_loss$Protein_category <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5, "Associated in both studies",
                                       ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y > 6.2e-5, "Associated in DiRECT only",
                                              ifelse(weight_loss$P.x > 2.1e-5 & weight_loss$P.y < 6.2e-5, "Associated in By-Band-Sleeve only", "Not associated")))

lm<- summary(lm(data=weight_loss, Fixed_effects_estimate ~ `Intervention*v2_Beta`))
#pdf(file = "DiRECT_BBS_comparison_pval_4groups_int_control_soma2.pdf", width = 12, height = 10)
ggplot(weight_loss, aes( x = `Intervention*v2_Beta`, y = Fixed_effects_estimate )) +
  geom_vline(xintercept=0, col="black", linetype = "dashed") +
  geom_hline(yintercept=0, col="black", linetype = "dashed")  +
  geom_point(aes(color= Protein_category), size = 1.5) +
  theme(legend.position = "bottom", legend.text = element_text(size=10)) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  geom_abline(intercept = 0, slope = 1, col="black", linetype = "dashed") +
  geom_errorbar(aes(ymin = Fixed_effects_estimate-SE,ymax = Fixed_effects_estimate+SE), alpha = 0.1) + 
  geom_errorbarh(aes(xmin = `Intervention*v2_Beta`-`Intervention*V2_SE`,xmax = `Intervention*v2_Beta`+`Intervention*V2_SE`), alpha = 0.1) +
  geom_label_repel(aes(label = ifelse( Protein_category == "Associated in both studies" , Target, "")), size = 4, color = "black", max.overlaps = 30) +
  xlab("Caloric restriction effect on protein levels in DiRECT (SDs)") + ylab("Bariatric surgery effect on protein levels in By-Band-Sleeve (SDs)") + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  annotate("text", x = -0.8, y = 1, label = "N=989", size = 5) +
  annotate("text", x = -0.8, y = 0.9, label = ("R squared = 0.05"), size = 5) +
  annotate("text", x = -0.8, y = 0.8, label = ("P = 1.1E-12"), size = 5) 
dev.off()

## Repeat the plot above with the correlation coefficients
correlation <- read_excel(path = "correlation_somalogic_olink.xlsx")
correlation <- distinct(correlation, SomaId, .keep_all = T)

results_correlation <- merge(weight_loss, correlation, by = "SomaId", all=F)
w <- which(abs(results_correlation$r) > 0.3)
results_correlation_plot <- results_correlation[w,]

lm<- summary(lm(data=results_correlation_plot, Fixed_effects_estimate ~ `Intervention*v2_Beta`))
#pdf(file = "DiRECT_BBS_comparison_pval_4groups_int_control_soma2_correlation.pdf", width = 12, height = 10)
ggplot(results_correlation_plot, aes( x = `Intervention*v2_Beta`, y = Fixed_effects_estimate )) +
  geom_vline(xintercept=0, col="black", linetype = "dashed") +
  geom_hline(yintercept=0, col="black", linetype = "dashed")  +
  geom_point(aes(color= Protein_category), size = 1.5) +
  theme(legend.position = "bottom", legend.text = element_text(size=10)) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  geom_abline(intercept = 0, slope = 1, col="black", linetype = "dashed") +
  geom_errorbar(aes(ymin = Fixed_effects_estimate-SE,ymax = Fixed_effects_estimate+SE), alpha = 0.1) + 
  geom_errorbarh(aes(xmin = `Intervention*v2_Beta`-`Intervention*V2_SE`,xmax = `Intervention*v2_Beta`+`Intervention*V2_SE`), alpha = 0.1) +
  geom_label_repel(aes(label = ifelse( Protein_category == "Associated in both studies" , Target, "")), size = 4, color = "black", max.overlaps = 20) +
  xlab("Caloric restriction effect on protein levels in DiRECT (SDs)") + ylab("Bariatric surgery effect on protein levels in By-Band-Sleeve (SDs)") + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"), legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  annotate("text", x = -0.8, y = 1, label = "N=423", size = 5) +
  annotate("text", x = -0.8, y = 0.9, label = ("R squared = 0.12"), size = 5) +
  annotate("text", x = -0.8, y = 0.8, label = ("P = 7.38E-12"), size = 5) 
dev.off()

# Comparison of estimates - only those that pass p value
w <- which(weight_loss$Group == "Associated in both studies")
weight_loss_consistent <- weight_loss[w,]
lm<- summary(lm(data=weight_loss_consistent, Fixed_effects_estimate ~ `Intervention*v2_Beta`))
#pdf(file = "DiRECT_BBS_comparison_consistent_pval.pdf", width = 12, height = 10)
ggplot(weight_loss_consistent, aes( x = `Intervention*v2_Beta`, y = Fixed_effects_estimate )) +
  geom_point() +
  geom_errorbar(aes(ymin = Fixed_effects_estimate-SE,ymax = Fixed_effects_estimate+SE), alpha = 0.1) + 
  geom_errorbarh(aes(xmin = `Intervention*v2_Beta`-`Intervention*V2_SE`,xmax = `Intervention*v2_Beta`+`Intervention*V2_SE`), alpha = 0.1) +
  geom_text(aes(label = Target), size = 3, color = "black") +
  geom_smooth(method = 'lm', formula = y ~ x) +
  xlab("DiRECT intervention effect on protein (SD)") + ylab("BBS intervention effect on protein (SD)") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) +
  annotate("text", x = -0.8, y = 1, label = "N=989", size = 5) +
  annotate("text", x = -0.8, y = 0.9, label = ("R squared = 0.89"), size = 5) +
  annotate("text", x = -0.8, y = 0.8, label = ("P = 3.2E-20"), size = 5)
dev.off()

### Merge in MR results to create forest plot of each of the 3 groups.
weight_loss <- merge(DiRECT_LMM, BBS_LMM, by="UniProt")
weight_loss$Association <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5,  "Associated in both studies", "Not associated in both studies")
weight_loss$Opposite_direction <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` < 0 & weight_loss$Fixed_effects_estimate > 0, "Opposite", 
                                         ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` > 0 & weight_loss$Fixed_effects_estimate < 0, "Opposite", "Not opposite"))
weight_loss$Group <- ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` < 0 & weight_loss$Fixed_effects_estimate < 0, "Associated in both studies",
                            ifelse(weight_loss$P.x < 2.1e-5 & weight_loss$P.y < 6.2e-5 & weight_loss$`Intervention*v2_Beta` > 0 & weight_loss$Fixed_effects_estimate > 0, "Associated in both studies", "Not associated in both studies"))
write.table(weight_loss, file = "DiRECT_BBS_results_soma2.txt", col.names = T, row.names = F, sep = "\t")

## Of those associated in the merged dataset, how many are negative and how many are positive
w <- which(weight_loss$P.x < 2.1e-5)
dir <- weight_loss[w,]
sum(dir$`Intervention*v2_Beta` < 0) ## 42
sum(dir$`Intervention*v2_Beta` > 0) ## 39

w <- which(weight_loss$P.y < 6.2e-5)
bbs <- weight_loss[w,]
sum(bbs$Fixed_effects_estimate < 0) ## 46
sum(bbs$Fixed_effects_estimate > 0) ## 84

w <- which(weight_loss$Association == "Not associated in both studies") ## 949 not associated in both studies

## Forest plot of the proteins with opposite directions, with MR results
w <- which(weight_loss$Opposite_direction == "Opposite")
Opp_data <- weight_loss[w,]
Opp_data <- merge(Opp_data, INTERVALMR, by.x = "UniProt", by.y="UniProtID", all.x=T, all.y=F)
Opp_data <- distinct(Opp_data,Target.x, .keep_all = T) 

## Format this data to make forest plot
Opp_forest <- as.data.frame(matrix(nrow=6, ncol=5))
colnames(Opp_forest) <- c("Protein", "Group", "Beta", "SE", "P")
Opp_forest$Group[1:2] <- "DiRECT"
Opp_forest$Group[3:4] <- "By-Band-Sleeve"
Opp_forest$Group[5:6]<- "INTERVAL"
Opp_forest$Group <- as.factor(Opp_forest$Group)
Opp_forest$Protein <- as.factor(Opp_forest$Protein)

Opp_forest$Protein <- Opp_data$Target.x
Opp_forest$Beta <- c(Opp_data$`Intervention*v2_Beta`, Opp_data$Fixed_effects_estimate, -Opp_data$Beta_coefficient)
Opp_forest$SE <- c(Opp_data$`Intervention*V2_SE`, Opp_data$SE.x, Opp_data$SE.y)
Opp_forest$P <- c(Opp_data$P.x, Opp_data$P.y, Opp_data$P_val)

#pdf(file = "forest_opposite_direction_soma2_v2.pdf", height = 8, width = 8)
forestplot(
  df = Opp_forest,
  estimate = Beta,
  se = SE,
  name = Protein,
  logodds = FALSE,
  pvalue= P,
  psignif =  0.05,
  ci=0.95,
  colour = Group,
  xlab = "Difference in protein (SDs ± 95% CI) per SD lower BMI (INTERVAL)\n or change in protein after intervention in DiRECT/By-Band-Sleeve (SDs ± 95% CI)"
)
dev.off()

## Proteins which increase after intervention in both studies
w <- which(weight_loss$Group == "Associated in both studies" &  weight_loss$`Intervention*v2_Beta` > 0 & weight_loss$Fixed_effects_estimate > 0)
Increased_proteins <- weight_loss[w,] 
Increased_proteins <- merge(Increased_proteins, INTERVALMR, by.x = "UniProt", by.y="UniProtID", all.x=T, all.y=F)
Increased_proteins <- distinct(Increased_proteins,Target.x, .keep_all = T) 

## Format this data to make forest plot
Increased_forest <- as.data.frame(matrix(nrow=39, ncol=5))
colnames(Increased_forest) <- c("Protein", "Group", "Beta", "SE", "P")
Increased_forest$Group[1:13] <- "DiRECT"
Increased_forest$Group[14:26] <- "By-Band-Sleeve"
Increased_forest$Group[27:39] <- "INTERVAL"
Increased_forest$Group <- as.factor(Increased_forest$Group)
Increased_forest$Protein <- as.factor(Increased_forest$Protein)

Increased_forest$Protein <- Increased_proteins$Target.x
Increased_forest$Beta <- c(Increased_proteins$`Intervention*v2_Beta`, Increased_proteins$Fixed_effects_estimate, -Increased_proteins$Beta_coefficient)
Increased_forest$SE <- c(Increased_proteins$`Intervention*V2_SE`, Increased_proteins$SE.x, Increased_proteins$SE.y)
Increased_forest$P <- c(Increased_proteins$P.x, Increased_proteins$P.y, Increased_proteins$P_val)

## Order by group then beta coefficient
Increased_forest$Group <- factor(Increased_forest$Group, levels = c("INTERVAL", "DiRECT", "By-Band-Sleeve"))
Increased_forest <- Increased_forest[order(Increased_forest$Group, -Increased_forest$Beta), ]

#pdf(file = "forest_increased_intervention_soma2_v3.pdf", height = 8, width = 8)
  p <- forestplot(
  df = Increased_forest,
  estimate = Beta,
  se = SE,
  name = Protein,
  logodds = FALSE,
  pvalue= P,
  psignif =  0.05,
  ci=0.95,
  colour = Group,
  xlab = "Difference in protein (SDs ± 95% CI) per SD lower BMI (INTERVAL)\n or change in protein after caloric restriction or bariatric surgery intervention (SDs ± 95% CI)"
)
  p + theme(axis.text = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))
dev.off()

## Proteins decreased with weight loss
w <- which(weight_loss$Group == "Associated in both studies" &  weight_loss$`Intervention*v2_Beta` < 0 & weight_loss$Fixed_effects_estimate < 0)
Decreased_proteins <- weight_loss[w,] 
Decreased_proteins <- merge(Decreased_proteins, INTERVALMR, by.x = "UniProt", by.y="UniProtID", all.x=T, all.y=F)
Decreased_proteins <- distinct(Decreased_proteins,Target.x, .keep_all = T) 

## Format this data to make forest plot
Decreased_forest <- as.data.frame(matrix(nrow=30, ncol=5))
colnames(Decreased_forest) <- c("Protein", "Group", "Beta", "SE", "P")
Decreased_forest$Group[1:10] <- "DiRECT"
Decreased_forest$Group[11:20] <- "By-Band-Sleeve"
Decreased_forest$Group[21:30] <- "INTERVAL"
Decreased_forest$Group <- as.factor(Decreased_forest$Group)
Decreased_forest$Protein <- as.factor(Decreased_forest$Protein)

Decreased_forest$Protein <- Decreased_proteins$Target.x
Decreased_forest$Beta <- c(Decreased_proteins$`Intervention*v2_Beta`, Decreased_proteins$Fixed_effects_estimate, -Decreased_proteins$Beta_coefficient)
Decreased_forest$SE <- c(Decreased_proteins$`Intervention*V2_SE`, Decreased_proteins$SE.x, Decreased_proteins$SE.y)
Decreased_forest$P <- c(Decreased_proteins$P.x, Decreased_proteins$P.y, Decreased_proteins$P_val)

Decreased_forest$Group <- factor(Decreased_forest$Group, levels = c("INTERVAL", "DiRECT", "By-Band-Sleeve"))
Decreased_forest <- Decreased_forest[order(Decreased_forest$Group, Decreased_forest$Beta), ]

#pdf(file = "forest_decreased_intervention_soma2_v3.pdf", height = 8, width = 8)
p <- forestplot(
  df = Decreased_forest,
  estimate = Beta,
  se = SE,
  name = Protein,
  logodds = FALSE,
  pvalue= P,
  psignif =  0.05,
  ci=0.95,
  colour = Group,
  xlab = "Difference in protein (SDs ± 95% CI) per SD lower BMI (INTERVAL)\n or change in protein after caloric restriction or bariatric surgery intervention (SDs ± 95% CI)"
)
p + theme(axis.text = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))
dev.off()

## Forest plot of just proteins with consistent evidence across three studies
x <- c("ADH4", "IGFBP1", "IGFBP2", "BCAN", "RTN4R", "IL1RN")
w <- which(weight_loss$Assay %in% x)
consistent_three <- weight_loss[w,]
consistent_three <- merge(consistent_three, INTERVALMR, by.x = "UniProt", by.y="UniProtID", all.x=T, all.y=F)
consistent_three <- distinct(consistent_three,Target.x, .keep_all = T) 

consistent_three_df <- as.data.frame(matrix(nrow=18, ncol=5))
colnames(consistent_three_df) <- c("Protein", "Group", "Beta", "SE", "P")
consistent_three_df$Group[1:6] <- "DiRECT"
consistent_three_df$Group[7:12] <- "By-Band-Sleeve"
consistent_three_df$Group[13:18] <- "INTERVAL"
consistent_three_df$Group <- as.factor(consistent_three_df$Group)
consistent_three_df$Protein <- as.factor(consistent_three_df$Protein)

consistent_three_df$Protein <- consistent_three$Target.x
consistent_three_df$Beta <- c(consistent_three$`Intervention*v2_Beta`, consistent_three$Fixed_effects_estimate, -consistent_three$Beta_coefficient)
consistent_three_df$SE <- c(consistent_three$`Intervention*V2_SE`, consistent_three$SE.x, consistent_three$SE.y)
consistent_three_df$P <- c(consistent_three$P.x, consistent_three$P.y, consistent_three$P_val)

#pdf(file = "forest_consistent_changed_3_soma2.pdf", height = 8, width = 8)
p <- forestplot(
  df = consistent_three_df,
  estimate = Beta,
  se = SE,
  name = Protein,
  logodds = FALSE,
  pvalue= P,
  psignif =  0.05,
  ci=0.95,
  colour = Group,
  xlab = "Difference in protein (SDs ± 95% CI) per SD lower BMI (INTERVAL)\n or change in protein after caloric restriction or bariatric surgery intervention (SDs ± 95% CI)"
)
p + theme(axis.text = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))
dev.off()

### Add in the correlation between somalogic and olink
correlation <- read_excel(path = "correlation_somalogic_olink.xlsx")
correlation <- distinct(correlation, SomaId, .keep_all = T)

results_correlation <- merge(weight_loss, correlation, by = "SomaId", all.x=T, all.y = F)
#write.xlsx(results_correlation, file = "Combined_results_with_correlations_soma2.xlsx", colNames=T, rowNames = F, sep = "\t")

## Mean and SD of correlation coefficient
w <- which(results_correlation$Group == "Associated in both studies")
corr_agreement <- results_correlation[w,]
mean(na.omit(corr_agreement$r))
min(na.omit(corr_agreement$r))
max(na.omit(corr_agreement$r))

## Plot of protein effect in DiRECT and BBS in proteins that are associated with both interventions
w <- which(results_correlation$Group == "Associated in both studies")
weight_loss_consistent_corr <- results_correlation[w,]
lm<- summary(lm(data=weight_loss_consistent_corr, Fixed_effects_estimate ~ `Intervention*v2_Beta`))
#pdf(file = "DiRECT_BBS_comparison_consistent_pval_soma2.pdf", width = 12, height = 10)
ggplot(weight_loss_consistent_corr, aes( x = `Intervention*v2_Beta`, y = Fixed_effects_estimate )) +
  geom_point(aes(color = r), size = 3) +
  geom_errorbar(aes(ymin = Fixed_effects_estimate-SE,ymax = Fixed_effects_estimate+SE), alpha = 0.1) + 
  geom_errorbarh(aes(xmin = `Intervention*v2_Beta`-`Intervention*V2_SE`,xmax = `Intervention*v2_Beta`+`Intervention*V2_SE`), alpha = 0.1) +
  geom_text(aes(label = Target), size = 3, color = "black") +
  geom_smooth(method = 'lm', formula = y ~ x) +
  xlab("DiRECT intervention effect on protein (SD)") + ylab("BBS intervention effect on protein (SD)") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) +
  annotate("text", x = -0.8, y = 1, label = "N=989", size = 5) +
  annotate("text", x = -0.8, y = 0.9, label = ("R squared = 0.89"), size = 5) +
  annotate("text", x = -0.8, y = 0.8, label = ("P = 3.2E-20"), size = 5)
dev.off()
