################################################################################

# Yasmijn Balder
# 23-06-2020

################################################################################

# Load all packages
library(missForest)

################################################################################

# Load the data
data <- read.csv("../Provided Data/Gene_cytokine_common.txt", 
                 sep = '\t')

# Select the annotation variables and make a new data frame
Annotation <- data[,c('Type_num', 'Sex_y', 'Comorbidity_y', "Amputation_y", 
                      "Septicshock_at_Baseline", "IVIG_y", "Clindamycin_y", 
                      'Death')]

colnames(Annotation)<-c('NSTI type','Sex','Comorbidity',"Clindamycin", "IVIG", 
                        "Septicshock at baseline", "Amputation", 'Death')

levels <- NULL
remove_levels <- NULL
# See how many levels each variable has
for (column in 1:ncol(Annotation)) {
  level <- length(unique(Annotation[,column]))
  if (level==1) {
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
Annotation <- Annotation[,-remove_levels]

# Deselect the annotation variables from the data
data <- subset(data, select = -c(ID,Type_num,Sex_y,Comorbidity_y,
                                 Septicshock_at_Baseline,IVIG_y,Clindamycin_y,
                                 Amputation_y,Death))

# Impute the missing data (default setting with 100 trees because small dataset)
data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE)
data_imp$OOBerror
data <- as.data.frame(data_imp$ximp)

# See how many levels each variable has
for (column in 1:ncol(data)) {
  level <- length(unique(data[,column]))
  if (level==1) {
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
data <- data[,-remove_levels]
# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################

# Perform t-tests

data <- cbind(data, Annotation)

significance <- NULL

################################################################################
# Clindamycin

clinda_0 <- data[which(data$Clindamycin==0),]
clinda_1 <- data[which(data$Clindamycin==1),]

# Compare ABCB5
p.value <- t.test(clinda_1$ABCB5, clinda_0$ABCB5, 
                  alternative = 'greater')$p.value
row <- c('Clindamycin', 'ABCB5', p.value)
significance <- cbind(significance, row)
# Compare EPHA3
p.value <- t.test(clinda_1$EPHA3, clinda_0$EPHA3, 
                  alternative = 'greater')$p.value
row <- c('Clindamycin', 'EPHA3', p.value)
significance <- cbind(significance, row)
# Compare LAMC3
p.value <- t.test(clinda_1$LAMC3, clinda_0$LAMC3, 
                  alternative = 'greater')$p.value
row <- c('Clindamycin', 'LAMC3', p.value)
significance <- cbind(significance, row)
# Compare SNAP91
p.value <- t.test(clinda_1$SNAP91, clinda_0$SNAP91, 
                  alternative = 'greater')$p.value
row <- c('Clindamycin', 'SNAP91', p.value)
significance <- cbind(significance, row)
# Compare Thrombomodulin
p.value <- t.test(clinda_1$Thrombomodulin, clinda_0$Thrombomodulin, 
                  alternative = 'greater')$p.value
row <- c('Clindamycin', 'Thrombomodulin', p.value)
significance <- cbind(significance, row)

################################################################################
# Comorbidity

como_0 <- data[which(data$Comorbidity==0),]
como_1 <- data[which(data$Comorbidity==1),]

# Compare PFKP
p.value <- t.test(como_1$PFKP, como_0$PFKP, 
                  alternative = 'greater')$p.value
row <- c('Comorbidity', 'PFKP', p.value)
significance <- cbind(significance, row)
# Compare ICAM1
p.value <- t.test(como_1$ICAM1, como_0$ICAM1, 
                  alternative = 'less')$p.value
row <- c('Comorbidity', 'ICAM1', p.value)
significance <- cbind(significance, row)

################################################################################
# Death

death_0 <- data[which(data$Death==0),]
death_1 <- data[which(data$Death==1),]

# Compare IL18
p.value <- t.test(death_1$IL18, death_0$IL18, alternative = 'greater')$p.value
row <- c('Death', 'IL18', p.value)
significance <- cbind(significance, row)
# Compare SELE
p.value <- t.test(death_1$SELE, death_0$SELE, 
                  alternative = 'greater')$p.value
row <- c('Death', 'SELE', p.value)
significance <- cbind(significance, row)
# Compare ARHGEF5
p.value <- t.test(death_1$ARHGEF5, death_0$ARHGEF5, 
                  alternative = 'greater')$p.value
row <- c('Death', 'ARHGEF5', p.value)
significance <- cbind(significance, row)
# Compare ACSM2B
p.value <- t.test(death_1$ACSM2B, death_0$ACSM2B, 
                  alternative = 'greater')$p.value
row <- c('Death', 'ACSM2B', p.value)
significance <- cbind(significance, row)

################################################################################
# IVIG

IVIG_0 <- data[which(data$IVIG==0),]
IVIG_1 <- data[which(data$IVIG==1),]

# Compare CXCL10IP10
p.value <- t.test(IVIG_1$CXCL10IP10, IVIG_0$CXCL10IP10, 
                  alternative = 'greater')$p.value
row <- c('IVIG', 'CXCL10IP10', p.value)
significance <- cbind(significance, row)
# Compare IL2
p.value <- t.test(IVIG_1$IL2, IVIG_0$IL2, alternative = 'greater')$p.value
row <- c('IVIG', 'IL2', p.value)
significance <- cbind(significance, row)
# Compare Ialpha1COL1A1
p.value <- t.test(IVIG_1$Ialpha1COL1A1, IVIG_0$Ialpha1COL1A1, 
                  alternative = 'greater')$p.value
row <- c('IVIG', 'Ialpha1COL1A1', p.value)
significance <- cbind(significance, row)
# Compare TENM1
p.value <- t.test(IVIG_1$TENM1, IVIG_0$TENM1, 
                  alternative = 'less')$p.value
row <- c('IVIG', 'TENM1', p.value)
significance <- cbind(significance, row)
# Compare CPS1
p.value <- t.test(IVIG_1$CPS1, IVIG_0$CPS1, 
                  alternative = 'less')$p.value
row <- c('IVIG', 'CPS1', p.value)
significance <- cbind(significance, row)
# Compare ADAMTS6
p.value <- t.test(IVIG_1$ADAMTS6, IVIG_0$ADAMTS6, 
                  alternative = 'less')$p.value
row <- c('IVIG', 'ADAMTS6', p.value)
significance <- cbind(significance, row)

################################################################################
# NSTI type

type_0 <- data[which(data$`NSTI type`==0),]
type_1 <- data[which(data$`NSTI type`==1),]

# Compare TG
p.value <- t.test(type_1$TG, type_0$TG, alternative = 'greater')$p.value
row <- c('NSTI type', 'TG', p.value)
significance <- cbind(significance, row)

################################################################################
# Septic shock

septicshock_0 <- data[which(data$Septicshock==0),]
septicshock_1 <- data[which(data$Septicshock==1),]

# Compare XK
p.value <- t.test(septicshock_1$XK, septicshock_0$XK, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'XK', p.value)
significance <- cbind(significance, row)
# Compare Resistin
p.value <- t.test(septicshock_1$Resistin, septicshock_0$Resistin, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'Resistin', p.value)
significance <- cbind(significance, row)
# Compare IL1ra
p.value <- t.test(septicshock_1$IL1ra, septicshock_0$IL1ra, 
                  alternative = 'greater')$p.value
row <- rbind.data.frame('Septic shock', 'IL1ra', p.value)
significance <- cbind(significance, row)
# Compare LGALS14
p.value <- t.test(septicshock_1$LGALS14, septicshock_0$LGALS14, 
                  alternative = 'less')$p.value
row <- rbind.data.frame('Septic shock', 'LGALS14', p.value)
significance <- cbind(significance, row)

################################################################################
# Sex

sex_0 <- data[which(data$Sex==0),]
sex_1 <- data[which(data$Sex==1),]

# Compare GABRA3
p.value <- t.test(sex_1$GABRA3, sex_0$GABRA3, alternative = 'less')$p.value
row <- c('Sex', 'GABRA3', p.value)
significance <- cbind(significance, row)
# Compare KDM5D
p.value <- t.test(sex_1$KDM5D, sex_0$KDM5D, alternative = 'greater')$p.value
row <- c('Sex', 'KDM5D', p.value)
significance <- cbind(significance, row)
# Compare NEXMIF
p.value <- t.test(sex_1$NEXMIF, sex_0$NEXMIF, alternative = 'less')$p.value
row <- c('Sex', 'NEXMIF', p.value)
significance <- cbind(significance, row)
# Compare DDX3Y
p.value <- t.test(sex_1$DDX3Y, sex_0$DDX3Y, alternative = 'greater')$p.value
row <- c('Sex', 'DDX3Y', p.value)
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
# write.csv(significance, 'P-values of genes+cytokines.csv')

# p1 <- hist(type_0$DeathOnly, breaks = 2)
# p2 <- hist(type_1$DeathOnly, breaks = 2)
# plot(p1, col=rgb(0,1,1,1/2))
# plot(p2, col=rgb(1,1,0,1/2), add=T)
# legend('topright',legend = c('Poly', 'Mono'), 
# col = c(rgb(0,1,1,1/2),rgb(1,1,0,1/2)), pch = 16)