################################################################################

# Load all packages
library(randomForest)
library(missForest)
library(gdata)
library(Rmisc)
library(svMisc)   # For progress bar, not strictly necessary

################################################################################
## LOAD THE DATA

# Load the data
data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Imputed_day0.csv")

################################################################################
## PREPARE THE DATA

# Remove some columns
remove<-NULL
for (column in c('IL1ra','IL23ELISA','IL33ELISA','ID_day','PatientID','day')){
  remove <- c(remove, which(colnames(data)==column))
}
data <- data[,-remove]

# Make a factor of the case type columns to be able to use it as a response
data$Case_type <- as.factor(data$Case_type)
levels(data$Case_type) <- c('Control','NSTI 1','NSTI 2','Cellulitis')
# Make a factor of the septic shock column
data$Septicshock <- as.factor(data$Septicshock)
levels(data$Septicshock) <- c('No', 'Yes')

# Which columns are discrete
discrete <- NULL
# Assign all factors to discrete
for (column in 1:ncol(data)) {
  if (class(data[,column])=='character') {
    discrete <- c(discrete, column)
  }
}
# Assign other discrete columns
for (column in c("Case_type","Sex","Comorbidity","Septicshock","IVIG",
                 "Clindamycin","AmputationOnly","DeathOnly","DeathAmputation",
                 "incl.amputation","inf.extremity")) {
  position <- which(colnames(data)==column)
  discrete <- c(discrete, position)
}
# The other columns are continuous
continuous <- c(1:ncol(data))
continuous <- continuous[-discrete]
# colnames(data)[continuous]

################################################################################
## PERFORM RANDOMFOREST

################################################################################
## Poly vs Mono (only NSTI samples)

# Select samples
data_NSTI <- data[which(data$Case_type!='Control'),]
data_NSTI <- data_NSTI[which(data_NSTI$Case_type!='Cellulitis'),]
# Make sure levels are dropped as well
data_NSTI <- drop.levels(data_NSTI)
# Keep only the response and the continuous variables
data_NSTI <- data_NSTI[,c(1, continuous)]

# Impute NAs in the continous columns
data_NSTI[,-1] <- missForest(xmis = as.matrix(data_NSTI[,-1]))$ximp

# See how many samples there are in each group
n1 = sum(data_NSTI$Case_type == 'NSTI 1')
n2 = sum(data_NSTI$Case_type == 'NSTI 2')

# Take 85% of the samples of the smallest group
nsize = round(min(n1,n2)*0.85) 

# Collect the error from the random forest 100 times
error_NSTI = NULL
for(k in 1 : 100){
  rf <- randomForest(data_NSTI$Case_type ~ ., data=data_NSTI, 
                     strata = data_NSTI$Case_type, sampsize = c(nsize,nsize))
  error_NSTI[k] = err.orig = mean(rf$err.rate[,1])
  progress(k)
}
CI_NSTI <- CI(error_NSTI,ci = 0.95)

################################################################################
# NSTI vs Controls

# Select samples
data_NSTI_control <- data[which(data$Case_type!='Cellulitis'),]
# Make one group of the NSTI types
levels(data_NSTI_control$Case_type) <- c('Control','NSTI','NSTI','Cellulitis')
# Make sure levels are dropped as well
data_NSTI_control <- drop.levels(data_NSTI_control)
# Keep only the response and the continuous variables
data_NSTI_control <- data_NSTI_control[,c(1, continuous)]

# Impute NAs in the continous columns
data_NSTI_control[,-1] <- missForest(xmis = as.matrix(data_NSTI_control[,-1]))$ximp

# See how many samples there are in each group
n1 = sum(data_NSTI_control$Case_type == 'NSTI')
n2 = sum(data_NSTI_control$Case_type == 'Control')

# Take 85% of the samples of the smallest group
nsize = round(min(n1,n2)*0.85) 

# Collect the error from the random forest 100 times
error_NSTI_control = NULL
for(k in 1 : 100){
  rf <- randomForest(data_NSTI_control$Case_type ~., data=data_NSTI_control, 
                     strata=data_NSTI_control$Case_type,sampsize=c(nsize,nsize))
  error_NSTI_control[k] = err.orig = mean(rf$err.rate[,1])
  progress(k)
}
CI_NSTI_control <- CI(error_NSTI_control,ci = 0.95)

################################################################################
# NSTI vs Non NSTI (case 3, cellulitis)

# Select samples
data_NSTI_non <- data[which(data$Case_type!='Control'),]
# Make one group of the NSTI types
levels(data_NSTI_non$Case_type) <- c('Control','NSTI','NSTI','Cellulitis')
# Make sure levels are dropped as well
data_NSTI_non <- drop.levels(data_NSTI_non)
# Keep only the response and the continuous variables
data_NSTI_non <- data_NSTI_non[,c(1, continuous)]

# Impute NAs in the continous columns
data_NSTI_non[,-1] <- missForest(xmis = as.matrix(data_NSTI_non[,-1]))$ximp

# See how many samples there are in each group
n1 = sum(data_NSTI_non$Case_type == 'NSTI')
n2 = sum(data_NSTI_non$Case_type == 'Cellulitis')

# Take 85% of the samples of the smallest group
nsize = round(min(n1,n2)*0.85) 

# Collect the error from the random forest 100 times
error_NSTI_non = NULL
for(k in 1 : 100){
  rf <- randomForest(data_NSTI_non$Case_type ~., data=data_NSTI_non, 
                     strata = data_NSTI_non$Case_type, sampsize=c(nsize,nsize))
  error_NSTI_non[k] = err.orig = mean(rf$err.rate[,1])
  progress(k)
}
CI_NSTI_non <- CI(error_NSTI_non,ci = 0.95)

################################################################################
# Non NSTI vs Controls

# Select samples
data_non_control <- data[which(data$Case_type!='NSTI 1'),]
data_non_control <- data_non_control[which(data_non_control$Case_type!='NSTI 2'),]
# Make sure levels are dropped as well
data_non_control <- drop.levels(data_non_control)
# Keep only the response and the continuous variables
data_non_control <- data_non_control[,c(1, continuous)]

# Impute NAs in the continous columns
data_non_control[,-1] <- as.data.frame(missForest(xmis=as.matrix(data_non_control[,-1]))$ximp)

# See how many samples there are in each group
n1 = sum(data_non_control$Case_type == 'Control')
n2 = sum(data_non_control$Case_type == 'Cellulitis')

# Take 85% of the samples of the smallest group
nsize = round(min(n1,n2)*0.85) 

# Collect the error from the random forest 100 times
error_non_control = NULL
for(k in 1 : 100){
  rf <- randomForest(data_non_control$Case_type ~ ., data=data_non_control, 
                     strata=data_non_control$Case_type, sampsize=c(nsize,nsize))
  error_non_control[k] = err.orig = mean(rf$err.rate[,1])
  progress(k)
}
CI_non_control <- CI(error_non_control,ci = 0.95)

################################################################################
# Septic shock vs non septic shock

# Select samples
data_ss <- data
# Keep only the response and the continuous variables
data_ss <- data_ss[,c(9, continuous)]

# Impute NAs in the continous columns
data_ss[,-1] <- missForest(xmis = as.matrix(data_ss[,-1]))$ximp
# Remove the NAs in the septic shock column
to_remove <- which(is.na(data_ss$Septicshock)==TRUE)
data_ss <- data_ss[-to_remove,]

# See how many samples there are in each group
n1 = sum(data_ss$Septicshock == 'No')
n2 = sum(data_ss$Septicshock == 'Yes')

# Take 85% of the samples of the smallest group
nsize = round(min(n1,n2)*0.85) 

# Collect the error from the random forest 100 times
error_ss = NULL
for(k in 1 : 100){
  rf <- randomForest(data_ss$Septicshock ~ ., data=data_ss, 
                     strata = data_ss$Septicshock, sampsize = c(nsize,nsize))
  error_ss[k] = err.orig = mean(rf$err.rate[,1])
  progress(k)
}
CI_ss <- CI(error_ss,ci = 0.95)

################################################################################
# Collect all errors

error_table <- rbind(CI_NSTI,CI_NSTI_control,CI_NSTI_non,CI_non_control,CI_ss)
rownames(error_table) <- c('Poly vs Mono','NSTI vs Controls','NSTI vs Non NSTI',
                      'Non NSTI vs Controls','Septic shock vs non septic shock')
View(error_table)

