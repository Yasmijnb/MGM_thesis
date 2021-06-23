###############################################################################

library(stringr)
setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")

###############################################################################

# Clinical (INFECT) + Gene data

clinical_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                          check.names = FALSE)
rownames(clinical_data) <- clinical_data$Row.names

gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/common_Gene.txt"
                      , sep = '\t', skip = 1)

rownames(gene_data) <- gene_data$ID
for (i in 1:nrow(gene_data)) {
  rownames(gene_data)[i] <- substring(rownames(gene_data)[i], 1, 4)
}

colnames(gene_data)[which(colnames(gene_data)=='Type_num')] <- 'NSTI type'
colnames(gene_data)[which(colnames(gene_data)=='Sex_y')] <- 'Sex'
colnames(gene_data)[which(colnames(gene_data)=='Comorbidity_y')] <- 'Comorbidity'
colnames(gene_data)[which(colnames(gene_data)=='Amputation_y')] <- 'Amputation'
colnames(gene_data)[which(colnames(gene_data)=='Septicshock_at_Baseline')] <- 'Septicshock'
colnames(gene_data)[which(colnames(gene_data)=='IVIG_y')] <- 'IVIG'
colnames(gene_data)[which(colnames(gene_data)=='Clindamycin_y')] <- 'Clindamycin'

data <- merge(clinical_data, gene_data, by.x = "row.names", by.y = "row.names")
rownames(data) <- data$'Row.names'
data <- data[,-c(1:2)]

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Clinical+Gene.csv')
data <- NULL

###############################################################################

# Clinical (INFECT) data + Gene & Cytokine data

clinical_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                          check.names = FALSE)
rownames(clinical_data) <- clinical_data$Row.names

gc_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Gene_cytokine_common.txt", 
                    sep = '\t')

rownames(gc_data) <- gc_data$ID
for (i in 1:nrow(gc_data)) {
  rownames(gc_data)[i] <- substring(rownames(gc_data)[i], 1,4)
}

colnames(gc_data)[which(colnames(gc_data)=='Type_num')] <- 'NSTI type'
colnames(gc_data)[which(colnames(gc_data)=='Sex_y')] <- 'Sex'
colnames(gc_data)[which(colnames(gc_data)=='Comorbidity_y')] <- 'Comorbidity'
colnames(gc_data)[which(colnames(gc_data)=='Amputation_y')] <- 'Amputation'
colnames(gc_data)[which(colnames(gc_data)=='Septicshock_at_Baseline')] <- 'Septicshock'
colnames(gc_data)[which(colnames(gc_data)=='IVIG_y')] <- 'IVIG'
colnames(gc_data)[which(colnames(gc_data)=='Clindamycin_y')] <- 'Clindamycin'

data <- merge(clinical_data, gc_data, by.x = "row.names", by.y = "row.names")
rownames(data) <- data$Row.names
data<-data[,-c(1:2)]

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Clinical+Gene+Cytokine.csv')
data <- NULL

###############################################################################

# Clinical (INFECT) + Cytokine data (OLD)

clinical_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                          check.names = FALSE)
rownames(clinical_data) <- clinical_data$Row.names

cytokines <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Final_Data_Imputed.csv")
cytokines <- subset(cytokines, select =-c(DO.NOT.USE..day,DO.NOT.USE.Case_type,
                                DO.NOT.USE.Microb_a,ï..ID_day))
cytokines <- cytokines[1:251,]
cytokines <- cytokines[-c(36,175),]
cytokines <- cytokines[,-which(colnames(cytokines)=="DeathAmputation")]
cytokines <- cytokines[,-which(colnames(cytokines)=="incl.amputation")]
cytokines <- cytokines[,-which(colnames(cytokines)=="inf.extremity")]
colnames(cytokines)[which(colnames(cytokines)=='Microbioogy')] <- 'NSTI type'
colnames(cytokines)[which(colnames(cytokines)=='DeathOnly')] <- 'Death'
colnames(cytokines)[which(colnames(cytokines)=='AmputationOnly')] <- 'Amputation'

data <- merge(clinical_data, cytokines, by.x = "row.names", by.y = "PatientID")
rownames(data) <- data$Row.names
data <- data[,-c(1:2)]

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Clinical+Cytokine(Old).csv')
data <- NULL

###############################################################################

# Clinical (INFECT) + Cytokine data (NEW)

clinical_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                          check.names = FALSE)
rownames(clinical_data) <- clinical_data$Row.names
clinical_data <- clinical_data[,-1]
clinical_data <- clinical_data[,-which(colnames(clinical_data)=='Row.names')]
clinical_data <- clinical_data[,-which(colnames(clinical_data)=='Patient ID')]
library("readxl")
cytokines <- read_excel("~/School/WUR/SSB-80336 - Thesis/Provided Data/Cytokines_3Dec2019.xlsx")
cytokines <- cytokines[,-which(colnames(cytokines)=="ID_Day")]
cytokines <- cytokines[,-which(colnames(cytokines)=="ID_day")]
cytokines$Case_type[which(cytokines$Case_type == 0)] <- 'Surgical control'
cytokines$Case_type[which(cytokines$Case_type == 1)] <- 'NSTI'
cytokines$Case_type[which(cytokines$Case_type == 2)] <- 'Non-NSTI'
cytokines$Case_type[which(cytokines$Case_type == 3)] <- 'Celullitis'
for (sample in 1:nrow(cytokines)) {
  if (cytokines$Case_type[sample]=='NSTI') {
    if (cytokines$Type_L[sample] == 'mono') {
      cytokines$Case_type[sample]<-'NSTI II'
    }
    if (cytokines$Type_L[sample] == 'poly') {
      cytokines$Case_type[sample]<-'NSTI I'
    }
  }
}
data <- merge(clinical_data, cytokines, by.x = "row.names", by.y = 'PatientID',
              all = TRUE, incomparables = NA, no.dups = FALSE)
data <- as.data.frame(data)
colnames(data)[which(colnames(data)=='Row.names')] <- 'PatientID'

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Clinical+Cytokine(New).csv')
data <- NULL

###############################################################################

# Bacterial and human genes

# Load the data
bact_gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/onlyBactFilt2ndProper.csv", 
                           sep = '\t')
human_gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/filt2ndHumanHgnc.csv", 
                            strip.white = TRUE, stringsAsFactors=FALSE, 
                            colClasses=c('character',rep('numeric', 102)))
clinic_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/metadataBiopsies_16s_Classified_ALL.csv")

# Change unknown into NA                                          ERROR
for (sample in 1:nrow(clinic_data)) {
  if (is.na(clinic_data$Tissue[sample])==FALSE) {
    if (clinic_data$Tissue[sample]=="unknown") {
      clinic_data$Tissue[sample] <- NA
    }
  }
  if (is.na(clinic_data$Tissue[sample])==FALSE) {
    if (clinic_data$Tissue[sample]=="") {
      clinic_data$Tissue[sample] <- NA
    }
  }
}

# Use the genes as row names
rownames(bact_gene_data) <- bact_gene_data[,1]
# Deselect the gene row
bact_gene_data <- bact_gene_data[,-1]
# Transpose the data
bact_gene_data <- as.data.frame(t(bact_gene_data))

# Use the genes as row names
rownames(human_gene_data) <- human_gene_data[,1]
# Deselect the gene row
human_gene_data <- human_gene_data[,-1]
# Transpose the data
human_gene_data <- as.data.frame(t(human_gene_data))

# Remove the X from the gene data so they match the clinical data
for (i in 1:nrow(human_gene_data)) {
  rownames(human_gene_data)[i] <- substring(rownames(human_gene_data)[i], 2)
}

# Remove the X from the gene data so they match the clinical data
for (i in 1:nrow(bact_gene_data)) {
  rownames(bact_gene_data)[i] <- substring(rownames(bact_gene_data)[i], 2)
}

# Change the colnames of the data into genes instead of UniRef90

library(stringr)  # For splitting the column name
library(UniprotR) # For retrieving the gene name from the uniref
library(svMisc)   # For progress bar, not strictly necessary

for (colnum in 1:ncol(bact_gene_data)) {
  # Parse the name
  splitted <- strsplit(colnames(bact_gene_data)[colnum], split = '_')[[1]]
  id_full <- splitted[2]
  id_num <- strsplit(id_full, split = '\\|')[[1]][1]
  # print(id_num)
  if (id_num != "unknown") {
    # If there is no match found, use the ID
    gene_name <- paste(id_num, '(UniprotID)')
    # Get the gene name
    tryCatch(gene_name <- GetProteinAnnontate(id_num, columns = 'entry name'),
             error = function(e) e)
    # Get the bacterium name
    bact <- paste(splitted[length(splitted)-1], splitted[length(splitted)])
    # Change the colname
    colnames(bact_gene_data)[colnum] <- paste0(gene_name,'; ', bact)
  }
  if (id_num == "unknown") {
    # Get the bacterium name
    bact <- paste(splitted[length(splitted)-1], splitted[length(splitted)])
    # Change the colname
    colnames(bact_gene_data)[colnum] <- paste("unknown;", bact)
  }
  # Keep track of progress
  progress(colnum*100/ncol(bact_gene_data),progress.bar = TRUE)
}

# Merge the data-frames
merge <- merge(bact_gene_data, clinic_data, by.x = "row.names", by.y = 'Sample',
               all = TRUE, incomparables = NA, no.dups = FALSE)
rownames(merge) <- merge$Row.names
merge <- merge[,-1]
data <- merge(human_gene_data, merge, by.x = "row.names", by.y = 'row.names',
              all = TRUE, incomparables = NA, no.dups = FALSE)
rownames(data) <- data$Row.names
data <- data[,-1]

# Select only day 1
data <- data[which(data$Day=='D1'),]

# Change column names for better understanding
colnames(data)[which(colnames(data) == 'Case')] <- 'City'
colnames(data)[which(colnames(data) == 'Predominant_species_humann2')] <- 
  'Predominant_species'
colnames(data)[which(colnames(data) == 'Microbioogy')] <- 'NSTI_type'
colnames(data)[which(colnames(data) == 'AmputationOnly')] <- 'Amputation'
colnames(data)[which(colnames(data) == 'DeathOnly')] <- 'Death'

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Bact+Human_genes.csv')
data <- NULL

###############################################################################

# All data

gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Bact+Human_genes.csv",
                      check.names = FALSE)
rownames(gene_data) <- gene_data[,1]
gene_data <- gene_data[,-1]
gene_data$SampleID <- rep(0, nrow(gene_data))
for (row in 1:nrow(gene_data)) {
  gene_data$SampleID[row] <- substring(rownames(gene_data)[row], 1, 4)
}
library("readxl")
cytokines <- read_excel("~/School/WUR/SSB-80336 - Thesis/Provided Data/Cytokines_3Dec2019.xlsx")
cytokines <- cytokines[,-which(colnames(cytokines)=="ID_Day")]
cytokines <- cytokines[,-which(colnames(cytokines)=="ID_day")]
cytokines$Case_type[which(cytokines$Case_type == 0)] <- 'Surgical control'
cytokines$Case_type[which(cytokines$Case_type == 1)] <- 'NSTI'
cytokines$Case_type[which(cytokines$Case_type == 2)] <- 'Non-NSTI'
cytokines$Case_type[which(cytokines$Case_type == 3)] <- 'Celullitis'
for (sample in 1:nrow(cytokines)) {
  if (cytokines$Case_type[sample]=='NSTI') {
    if (cytokines$Type_L[sample] == 'mono') {
      cytokines$Case_type[sample]<-'NSTI II'
    }
    if (cytokines$Type_L[sample] == 'poly') {
      cytokines$Case_type[sample]<-'NSTI I'
    }
  }
}
colnames(cytokines)[which(colnames(cytokines) == 'DeathOnly')] <- 'Death'
colnames(cytokines)[which(colnames(cytokines) == 'AmputationOnly')] <- 'Amputation'
cytokines <- cytokines[cytokines$day == 0,]
cytokines <- cytokines[cytokines$Case_type == 'NSTI I' | cytokines$Case_type == 'NSTI II',]
for (column in c('Microbiology','Type_L','Microb_a','SOFA2','SOFA3','SOFA4',
                 'SOFA5','SOFA6','SOFA7','incl.amputation','DeathAmputation','day')) {
  cytokines <- cytokines[,-which(colnames(cytokines)==column)]
}

# cytokine_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Final_Data_Imputed.csv")
# cytokine_data <- subset(cytokine_data, select = -c(DO.NOT.USE..day,
#                                                    DO.NOT.USE.Case_type,
#                                                    DO.NOT.USE.Microb_a,
#                                                    ï..ID_day,DeathAmputation,
#                                                    incl.amputation,
#                                                    inf.extremity))
# cytokine_data <- cytokine_data[1:251,]
# cytokine_data <- cytokine_data[-c(36,175),]
# colnames(cytokine_data)[which(colnames(cytokine_data) == 'Microbioogy')] <- 'NSTI_type'

# Merge
data <- merge(cytokines, gene_data, by.x = "PatientID", by.y = 'SampleID',
              all = TRUE, incomparables = NA, no.dups = FALSE)
data <- as.data.frame(data)

# When there is no entry, change it to NA
for (column in 1:ncol(data)) {
  levels(data[,column])[levels(data[,column]) == ""]<-NA
  # New R version:
  empty <- which(data[, column]=="")
  data[empty, column] <- NA
}

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'All_data.csv')
data <- NULL

###############################################################################

# All + INFECT

all_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/All_data.csv",
                     check.names = FALSE)
all_data <- all_data[,-1]
clinical_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                          check.names = FALSE)
rownames(clinical_data) <- clinical_data[,1]
clinical_data <- clinical_data[,-1]

data <- merge(all_data, clinical_data, by.x = "PatientID", by.y = 'row.names',
              all = TRUE, incomparables = NA, no.dups = FALSE)
data <- as.data.frame(data)

# When there is no entry, change it to NA
for (column in 1:ncol(data)) {
  levels(data[,column])[levels(data[,column]) == ""]<-NA
  # New R version:
  empty <- which(data[, column]=="")
  data[empty, column] <- NA
}

# Chance some names and remove some duplicates
data <- data[,-which(colnames(data)=='Row.names')]
data <- data[,-which(colnames(data)=='PatientId')]
data <- data[,-which(colnames(data)=='PatientID')]
data <- data[,-which(colnames(data)=='BMI.x')]
colnames(data)[which(colnames(data)=='BMI.y')] <- 'BMI'
data <- data[,-which(colnames(data)=='Age')]
colnames(data)[which(colnames(data)=='Death')] <- 'Death (0=no, 1=yes)'
data <- data[,-which(colnames(data)=='Amputation')]
data <- data[,-which(colnames(data)=='Septicshock')]
data <- data[,-which(colnames(data)=='Comorbidity')]
data <- data[,-which(colnames(data)=='City')]
data$`City (2=Copenhagen, 3=Stockholm, 4=Karlskrona, 5=Gothenburg, 6=Bergen)`[401] <- 6
colnames(data)[which(colnames(data)=='IVIG')] <- 'IVIG (0=no, 1=yes)'
data <- data[,-which(colnames(data)=='Sex')]
data <- data[,-which(colnames(data)=='Patient Dead Day 90 (0=no, 1=yes).1')]

setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'All_data+INFECT.csv')
data <- NULL

