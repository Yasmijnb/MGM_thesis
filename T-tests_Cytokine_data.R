################################################################################

# Load all packages
library(missForest)

################################################################################

# Load the data
data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Final_Data_Imputed.csv")

# Remove DO.NOT.USE columns
data <- subset(data, select =-c(DO.NOT.USE..day,DO.NOT.USE.Case_type,
                                DO.NOT.USE.Microb_a,ï..ID_day,PatientID))
# Remove samples without clinical data
data <- data[1:251,]      # The last few samples have no clinical data
data <- data[-c(36,175),] # These two samples have NA in discrete columns

# Select the annotation variables and make a new data frame
Annotation <- data[,c('Microbioogy','Sex','Comorbidity',"Septicshock","IVIG",
                      "Clindamycin","AmputationOnly",'DeathOnly',
                      "incl.amputation","inf.extremity")]

# Deselect the annotation variables from the data
data <- subset(data, select = -c(Microbioogy,Sex,Comorbidity,Septicshock, 
                                 IVIG,Clindamycin,AmputationOnly,DeathOnly, 
                                 DeathAmputation,incl.amputation,inf.extremity))

# Impute the missing data (default setting with 100 trees because small dataset)
data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE)
data_imp$OOBerror
data <- as.data.frame(data_imp$ximp)

# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################
# Perform t-tests

data <- cbind(data, Annotation)

significance <- NULL

################################################################################
# Amputation

amputation_0 <- data[which(data$AmputationOnly==0),]
amputation_1 <- data[which(data$AmputationOnly==1),]

# Compare C5C5a
p.value <- t.test(amputation_1$C5C5a, amputation_0$C5C5a, 
                  alternative = 'less')$p.value
row <- c('Amputation', 'C5C5a', p.value)
significance <- cbind(significance, row)
# Compare CollagenIValpha1
p.value <- t.test(amputation_1$CollagenIValpha1, amputation_0$CollagenIValpha1, 
                  alternative = 'greater')$p.value
row <- c('Amputation', 'CollagenIValpha1', p.value)
significance <- cbind(significance, row)
# Compare IL17A
p.value <- t.test(amputation_1$IL17A, amputation_0$IL17A, 
                  alternative = 'greater')$p.value
row <- c('Amputation', 'IL17A', p.value)
significance <- cbind(significance, row)
# Compare Thrombomodulin
p.value <- t.test(amputation_1$Thrombomodulin, amputation_0$Thrombomodulin, 
                  alternative = 'greater')$p.value
row <- c('Amputation', 'Thrombomodulin', p.value)
significance <- cbind(significance, row)
# Compare VCAM1
p.value <- t.test(amputation_1$VCAM1, amputation_0$VCAM1, 
                  alternative = 'greater')$p.value
row <- c('Amputation', 'VCAM1', p.value)
significance <- cbind(significance, row)
# Compare CCL5RANTES
p.value <- t.test(amputation_1$CCL5RANTES, amputation_0$CCL5RANTES, 
                  alternative = 'less')$p.value
row <- c('Amputation', 'CCL5RANTES', p.value)
significance <- cbind(significance, row)

################################################################################
# Comorbidity

comorbidity_0 <- data[which(data$Comorbidity==0),]
comorbidity_1 <- data[which(data$Comorbidity==1),]

# Compare ESelectin
p.value <- t.test(comorbidity_1$ESelectin, comorbidity_0$ESelectin, 
                  alternative = 'less')$p.value
row <- c('Comorbidity', 'ESelectin', p.value)
significance <- cbind(significance, row)
# Compare ICAM1
p.value <- t.test(comorbidity_1$ICAM1, comorbidity_0$ICAM1, 
                  alternative = 'less')$p.value
row <- c('Comorbidity', 'ICAM1', p.value)
significance <- cbind(significance, row)
# Compare S100A8
p.value <- t.test(comorbidity_1$S100A8, comorbidity_0$S100A8, 
                  alternative = 'less')$p.value
row <- c('Comorbidity', 'S100A8', p.value)
significance <- cbind(significance, row)

################################################################################
# Death

death_0 <- data[which(data$DeathOnly==0),]
death_1 <- data[which(data$DeathOnly==1),]

# Compare C5C5a
p.value <- t.test(death_1$C5C5a, death_0$C5C5a, alternative = 'less')$p.value
row <- c('Death', 'C5C5a', p.value)
significance <- cbind(significance, row)
# Compare CXCL10IP10
p.value <- t.test(death_1$CXCL10IP10, death_0$CXCL10IP10, 
                  alternative = 'less')$p.value
row <- c('Death', 'CXCL10IP10', p.value)
significance <- cbind(significance, row)
# Compare MMP1
p.value <- t.test(death_1$MMP1, death_0$MMP1, alternative = 'greater')$p.value
row <- c('Death', 'MMP1', p.value)
significance <- cbind(significance, row)
# Compare CCL5RANTES
p.value <- t.test(death_1$CCL5RANTES, death_0$CCL5RANTES, 
                  alternative = 'less')$p.value
row <- c('Death', 'CCL5RANTES', p.value)
significance <- cbind(significance, row)
# Compare Ialpha1COL1A1
p.value <- t.test(death_1$Ialpha1COL1A1, death_0$Ialpha1COL1A1, 
                  alternative = 'less')$p.value
row <- c('Death', 'Ialpha1COL1A1', p.value)
significance <- cbind(significance, row)

################################################################################
# IVIG

IVIG_0 <- data[which(data$IVIG==0),]
IVIG_1 <- data[which(data$IVIG==1),]

# Compare Resistin
p.value <- t.test(IVIG_1$Resistin, IVIG_0$Resistin, 
                  alternative = 'greater')$p.value
row <- c('IVIG', 'Resistin', p.value)
significance <- cbind(significance, row)
# Compare ICAM1
p.value <- t.test(IVIG_1$ICAM1, IVIG_0$ICAM1, alternative = 'greater')$p.value
row <- c('IVIG', 'ICAM1', p.value)
significance <- cbind(significance, row)
# Compare IL1ra
p.value <- t.test(IVIG_1$IL1ra, IVIG_0$IL1ra, alternative = 'greater')$p.value
row <- c('IVIG', 'IL1ra', p.value)
significance <- cbind(significance, row)
# Compare Thrombomodulin
p.value <- t.test(IVIG_1$Thrombomodulin, IVIG_0$Thrombomodulin, 
                  alternative = 'less')$p.value
row <- c('IVIG', 'Thrombomodulin', p.value)
significance <- cbind(significance, row)
# Compare IL23ELISA
p.value <- t.test(IVIG_1$IL23ELISA, IVIG_0$IL23ELISA, 
                  alternative = 'greater')$p.value
row <- c('IVIG', 'IL23ELISA', p.value)
significance <- cbind(significance, row)

################################################################################
# NSTI type

type_0 <- data[which(data$Microbioogy=='poly'),]
type_1 <- data[which(data$Microbioogy=='mono'),]

# Compare CXCL10IP10
p.value <- t.test(type_1$CXCL10IP10, type_0$CXCL10IP10, 
                  alternative = 'greater')$p.value
row <- c('NSTI type', 'CXCL10IP10', p.value)
significance <- cbind(significance, row)
# Compare ESelectin
p.value <- t.test(type_1$ESelectin, type_0$ESelectin, 
                  alternative = 'greater')$p.value
row <- c('NSTI type', 'ESelectin', p.value)
significance <- cbind(significance, row)
# Compare GCSF
p.value <- t.test(type_1$GCSF, type_0$GCSF, alternative = 'greater')$p.value
row <- c('NSTI type', 'GCSF', p.value)
significance <- cbind(significance, row)
# Compare IL2
p.value <- t.test(type_1$IL2, type_0$IL2, alternative = 'greater')$p.value
row <- c('NSTI type', 'IL2', p.value)
significance <- cbind(significance, row)
# Compare MMP8
p.value <- t.test(type_1$MMP8, type_0$MMP8, alternative = 'less')$p.value
row <- c('NSTI type', 'MMP8', p.value)
significance <- cbind(significance, row)
# Compare CCL2MCP1
p.value <- t.test(type_1$CCL2MCP1, type_0$CCL2MCP1, 
                  alternative = 'less')$p.value
row <- c('NSTI type', 'CCL2MCP1', p.value)
significance <- cbind(significance, row)
# Compare FasLigand
p.value <- t.test(type_1$FasLigand, type_0$FasLigand, 
                  alternative = 'greater')$p.value
row <- c('NSTI type', 'FasLigand', p.value)
significance <- cbind(significance, row)
# Compare S100A8
p.value <- t.test(type_1$S100A8, type_0$S100A8, alternative = 'greater')$p.value
row <- c('NSTI type', 'S100A8', p.value)
significance <- cbind(significance, row)
# Compare VCAM1
p.value <- t.test(type_1$VCAM1, type_0$VCAM1, alternative = 'greater')$p.value
row <- c('NSTI type', 'VCAM1', p.value)
significance <- cbind(significance, row)
# Compare MMP9
p.value <- t.test(type_1$MMP9, type_0$MMP9, alternative = 'less')$p.value
row <- c('NSTI type', 'MMP9', p.value)
significance <- cbind(significance, row)

################################################################################
# Septic shock

septicshock_0 <- data[which(data$Septicshock==0),]
septicshock_1 <- data[which(data$Septicshock==1),]

# Compare C5C5a
p.value <- t.test(septicshock_1$C5C5a, septicshock_0$C5C5a, 
                  alternative = 'less')$p.value
row <- rbind.data.frame('Septic shock', 'C5C5a', p.value)
significance <- cbind(significance, row)
# Compare IL2
p.value <- t.test(septicshock_1$IL2, septicshock_0$IL2, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'IL2', p.value)
significance <- cbind(significance, row)
# Compare IL36betaIL1F8
p.value <- t.test(septicshock_1$IL36betaIL1F8, septicshock_0$IL36betaIL1F8, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'IL36betaIL1F8', p.value)
significance <- cbind(significance, row)
# Compare CollagenIValpha1
p.value <- t.test(septicshock_1$CollagenIValpha1,septicshock_0$CollagenIValpha1, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'CollagenIValpha1', p.value)
significance <- cbind(significance, row)
# Compare IL1ra
p.value <- t.test(septicshock_1$IL1ra, septicshock_0$IL1ra, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'IL1ra', p.value)
significance <- cbind(significance, row)
# Compare IL4
p.value <- t.test(septicshock_1$IL4, septicshock_0$IL4, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'IL4', p.value)
significance <- cbind(significance, row)

################################################################################
# Sex

sex_0 <- data[which(data$Sex==0),]
sex_1 <- data[which(data$Sex==1),]

# Compare VCAM1
p.value <- t.test(sex_1$VCAM1, sex_0$VCAM1, alternative = 'greater')$p.value
row <- c('Sex', 'VCAM1', p.value)
significance <- cbind(significance, row)

################################################################################
# Significance 

# Transpose the data frame
significance <- t(significance)
# Add column names
colnames(significance) <- c("Clinical parameter","Cytokine","P-value")
# Get rid of the rownames
rownames(significance) <- 1:nrow(significance)

# Make the table into a dataframe
significance <- as.data.frame(significance)
# Make the p-value column integers
significance$`P-value` <- as.numeric(as.character(significance$`P-value`))

# Do multiple tesing correction
significance$BH <- p.adjust(significance$`P-value`, method = "BH")

# Make a conclusion column with 'significant'
significance$Conclusion <- rep('significant',nrow(significance))
# If the p-value is too high, change the conclusion to 'not significant'
for (ttest in 1:nrow(significance)) {
  if (significance[ttest,'BH'] > 0.05) {
    significance[ttest,'Conclusion'] <- 'not significant'
  }
}

# Round the P-values
significance$`P-value` <- round(significance$`P-value`, 5)
significance$BH <- round(significance$BH, 5)

View(significance)

# Save the table
# myPath <- "~/School/WUR/SSB-80336 - Thesis/"
# setwd(myPath)
# write.csv(significance, 'P-values of cytokines.csv')