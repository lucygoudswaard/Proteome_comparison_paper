##### Preparation of data for metaboprep
library(readxl)
library(tidyr)
library(psych)
library(dplyr)

## Set working directory
# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## Read in data
mydata <- read_excel(paste0(prot_release, "2022-05-19/data/raw/20202265_Timpson_NPX_2021-04-27_V2.xlsx"))

### Expore the data - how may assay warnings? Will not remove them for now.
w <- which(mydata$Assay_Warning == "WARN")
length(w) # 90 assay warnings
warn <- mydata[w,]

## How many QC warnings?
w <- which(mydata$QC_Warning == "WARN")
length(w) # 28422 assay warnings
qcwarn <- mydata[w,]

## How many proteins
length(unique(mydata$UniProt)) ## 1463
length(unique(mydata$SampleID)) ## 258
length(unique(mydata$OlinkID)) ## 1472

### Write data to run through metaboprep - restrict to sample ID and protein NPX
protein_data <- pivot_wider(
  mydata,
  id_cols = SampleID,
  names_from = OlinkID,
  names_prefix = "",
  names_sep = "",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = NPX,
  values_fill = NA,
  values_fn = NULL
)

#### Ensure data is in a dataframe
protein_data <- as.data.frame(protein_data)

#### First column names needs to be "id"
colnames(protein_data)[1] <- "id"

#### Edit sampleIDs so they do not have characters or uppercase letters
protein_data[,1] = gsub("_", "", protein_data[,1])
protein_data[,1] = gsub("-", "", protein_data[,1])
protein_data[,1] = tolower(protein_data[,1])

#### Change sampleID and olink IDs to be lower case
colnames(protein_data)[1:ncol(protein_data)] = tolower(colnames(protein_data)[1:ncol(protein_data)]) 
#write.table(protein_data, file = paste0(data_intermediate_dir, "olink_metaboprep.txt"), sep = "\t", col.names = T, row.names = F)
