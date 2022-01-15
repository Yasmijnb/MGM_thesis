################################################################################

# Yasmijn Balder
# 06-06-2020

################################################################################

# Load all packages
library(missForest)

################################################################################

# Load the data
data <- read.csv("../Provided Data/common_Gene.txt", sep = '\t', skip = 1)

# Select the annotation variables and make a new data frame
Annotation <- data[,c('Type_num','Sex_y','Comorbidity_y',"Clindamycin_y",
                      "Septicshock_at_Baseline","IVIG_y","Amputation_y",
                      'Death')]
                      

colnames(Annotation)<-c('NSTI type','Sex','Comorbidity',"Clindamycin",
                        "Septicshock at baseline","IVIG","Amputation",'Death')

levels <- NULL
remove_levels <- NULL
# See how many levels each variable has
for (column in 1:ncol(Annotation)) {
  level <- length(unique(Annotation[,column]))
  levels <- c(levels, level)
  if (level==1) {
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
Annotation <- Annotation[,-remove_levels]
levels <- levels[-remove_levels]

# Deselect the annotation variables from the data
data <- subset(data, select = -c(ID,Type_num,Sex_y,Comorbidity_y,Amputation_y,
                                 Septicshock_at_Baseline,IVIG_y,Clindamycin_y,
                                 Death))

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

amputation_0 <- data[which(data$Amputation==0),]
amputation_1 <- data[which(data$Amputation==1),]

# Compare ABCB5
p.value <- t.test(amputation_1$ABCB5, amputation_0$ABCB5, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Amputation', 'ABCB5', p.value)
significance <- row # Start like this, otherwise errors occur

# Compare EPHA3
p.value <- t.test(amputation_1$EPHA3, amputation_0$EPHA3, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Amputation', 'EPHA3', p.value)
significance <- cbind(significance, row)

# Compare LAMC3
p.value <- t.test(amputation_1$LAMC3, amputation_0$LAMC3, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Amputation', 'LAMC3', p.value)
significance <- cbind(significance, row)

# Compare COL23A1
p.value <- t.test(amputation_1$COL23A1, amputation_0$COL23A1, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Amputation', 'COL23A1', p.value)
significance <- cbind(significance, row)

# Compare SNAP91
p.value <- t.test(amputation_1$SNAP91, amputation_0$SNAP91, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Amputation', 'SNAP91', p.value)
significance <- cbind(significance, row)

################################################################################
# Death

death_0 <- data[which(data$Death==0),]
death_1 <- data[which(data$Death==1),]

# Compare PON1
p.value <- t.test(death_1$PON1, death_0$PON1, alternative = 'greater')$p.value
row <- rbind.data.frame('Death', 'PON1', p.value)
significance <- cbind(significance, row)

# Compare SELE
p.value <- t.test(death_1$SELE, death_0$SELE, alternative = 'greater')$p.value
row <- rbind.data.frame('Death', 'SELE', p.value)
significance <- cbind(significance, row)

# Compare ARHGEF5
p.value <- t.test(death_1$ARHGEF5, death_0$ARHGEF5, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Death', 'ARHGEF5', p.value)
significance <- cbind(significance, row)

# Compare ACSM2B
p.value <- t.test(death_1$ACSM2B, death_0$ACSM2B, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Death', 'ACSM2B', p.value)
significance <- cbind(significance, row)

################################################################################
# IVIG

IVIG_0 <- data[which(data$IVIG==0),]
IVIG_1 <- data[which(data$IVIG==1),]

# Compare LGALS14
p.value <- t.test(IVIG_1$LGALS14, IVIG_0$LGALS14, alternative = 'less')$p.value
row <- rbind.data.frame('IVIG', 'LGALS14', p.value)
significance <- cbind(significance, row)

# Compare XK
p.value <- t.test(IVIG_1$XK, IVIG_0$XK, alternative = 'greater')$p.value
row <- rbind.data.frame('IVIG', 'XK', p.value)
significance <- cbind(significance, row)

# Compare TMEM159
p.value <- t.test(IVIG_1$TMEM159, IVIG_0$TMEM159, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('IVIG', 'TMEM159', p.value)
significance <- cbind(significance, row)

# Compare SLC4A8
p.value <- t.test(IVIG_1$SLC4A8, IVIG_0$SLC4A8, alternative = 'less')$p.value
row <- rbind.data.frame('IVIG', 'SLC4A8', p.value)
significance <- cbind(significance, row)

################################################################################
# NSTI Type

NSTIt_0 <- data[which(data$`NSTI type`==0),]
NSTIt_1 <- data[which(data$`NSTI type`==1),]

# Compare TG
p.value <- t.test(NSTIt_1$TG, NSTIt_0$TG, alternative = 'greater')$p.value
row <- rbind.data.frame('NSTI Type', 'TG', p.value)
significance <- cbind(significance, row)

################################################################################
# Septic shock

septicshock_0 <- data[which(data$`Septicshock at baseline`==0),]
septicshock_1 <- data[which(data$`Septicshock at baseline`==1),]

# Compare TENM1
p.value <- t.test(septicshock_1$TENM1, septicshock_0$TENM1, 
                  alternative = 'less')$p.value
row <- rbind.data.frame('Septic shock', 'TENM1', p.value)
significance <- cbind(significance, row)

# Compare CPS1
p.value <- t.test(septicshock_1$CPS1, septicshock_0$CPS1, 
                  alternative = 'less')$p.value
row <- c('Septic shock', 'CPS1', p.value)
significance <- cbind(significance, row)

# Compare TG
p.value <- t.test(septicshock_1$TG, septicshock_0$TG, 
                  alternative = 'less')$p.value
row <- rbind.data.frame('Septic shock', 'TG', p.value)
significance <- cbind(significance, row)

# Compare ADAMTS6
p.value <- t.test(septicshock_1$ADAMTS6, septicshock_0$ADAMTS6, 
                  alternative = 'less')$p.value
row <- rbind.data.frame('Septic shock', 'ADAMTS6', p.value)
significance <- cbind(significance, row)

################################################################################
# Sex

sex_0 <- data[which(data$Sex==0),]
sex_1 <- data[which(data$Sex==1),]

# Compare GABRA3
p.value <- t.test(sex_1$GABRA3, sex_0$GABRA3, alternative = 'less')$p.value
row <- rbind.data.frame('Sex', 'GABRA3', p.value)
significance <- cbind(significance, row)

# Compare KDM5D
p.value <- t.test(sex_1$KDM5D, sex_0$KDM5D, alternative = 'greater')$p.value
row <- rbind.data.frame('Sex', 'KDM5D', p.value)
significance <- cbind(significance, row)

# Compare NEXMIF
p.value <- t.test(sex_1$NEXMIF, sex_0$NEXMIF, alternative = 'less')$p.value
row <- rbind.data.frame('Sex', 'NEXMIF', p.value)
significance <- cbind(significance, row)

# Compare DDX3Y
p.value <- t.test(sex_1$DDX3Y, sex_0$DDX3Y, alternative = 'greater')$p.value
row <- rbind.data.frame('Sex', 'DDX3Y', p.value)
significance <- cbind(significance, row)

################################################################################
# Significance 
# Transpose the data frame
significance <- t(significance)
# Add column names
colnames(significance) <- c("Clinical parameter","Gene","P-value")
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
# write.csv(significance, '../P-values of genes.csv')
