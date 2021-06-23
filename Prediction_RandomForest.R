################################################################################

# Yasmijn Balder
# 04-09-2020

# Predict with random forest

###############################################################################

# Load all packages
library(randomForest)
library(Rmisc)    # For CI (confidence interval)
library(svMisc)   # For progress, not strictly necessary
library(ggplot2)  # For the plot with mean decease gini

################################################################################

# Load the data
library("readxl")
data <- read_excel("~/School/WUR/SSB-80336 - Thesis/Provided Data/Cytokines_3Dec2019.xlsx")
data <- data[,-which(colnames(data)=="ID_Day")]
data <- data[,-which(colnames(data)=="ID_day")]
data <- data[,-which(colnames(data)=="day")]
data$Case_type[which(data$Case_type == 0)] <- 'Surgical control'
data$Case_type[which(data$Case_type == 1)] <- 'NSTI'
data$Case_type[which(data$Case_type == 2)] <- 'Non-NSTI'
data$Case_type[which(data$Case_type == 3)] <- 'Celullitis'
for (sample in 1:nrow(data)) {
  if (data$Case_type[sample]=='NSTI') {
    if (data$Type_L[sample] == 'mono') {
      data$Case_type[sample]<-'NSTI II'
    }
    if (data$Type_L[sample] == 'poly') {
      data$Case_type[sample]<-'NSTI I'
    }
  }
}

################################################################################

# Keep only columns which have less than a limit missing data
remove <- NULL
limit <- 0 # percentage of missing data that is allowed per variable
for (i in 26:ncol(data)){
  percentage <- 100*sum(is.na(data[,i]))/nrow(data)
  if (percentage>limit) {
    # print(percentage)
    remove <- c(remove, i)
  }
}
if (limit < 100) {
  data <- data[,-remove]
}

# Remove some columns
for (column in c('DeathAmputation','PatientID','SOFA2','SOFA3','SOFA4','SOFA5',
                 'SOFA6','SOFA7','incl.amputation','Microbiology','Type_L',
                 'Microb_a')) {
  data <- data[,-which(colnames(data)==column)]
}
# Change some column names
colnames(data)[which(colnames(data)=='DeathOnly')] <- 'Death'
colnames(data)[which(colnames(data)=='AmputationOnly')] <- 'Amputation'

################################################################################

# Patient charactaristics
n.male <- length(which(data$Sex == 1))
p.male <- n.male*100/nrow(data)
n.ss <- length(which(data$Septicshock == 1))
p.ss <- n.ss*100/nrow(data)
n.amp <- length(which(data$Amputation == 1))
p.amp <- n.amp*100/nrow(data)
n.death <- length(which(data$Death == 1))
p.death <- n.death*100/nrow(data)
n.ii <- length(which(data$Case_type=='NSTI II'))
p.ii <- n.ii*100/nrow(data)
n.i <- length(which(data$Case_type=='NSTI I'))
p.i <- n.i*100/nrow(data)
n.non <- length(which(data$Case_type=='Non-NSTI'))
p.non <- n.non*100/nrow(data)
n.cel <- length(which(data$Case_type=='Celullitis'))
p.cel <- n.cel*100/nrow(data)
n.control <- length(which(data$Case_type=='Surgical control'))
p.control <- n.cel*100/nrow(data)

pat.char <- rbind(c(n.amp, round(p.amp,2)),
                  c(n.death, round(p.death,2)),
                  c(n.ss, round(p.ss,2)),
                  c(n.male, round(p.male,2)),
                  c(n.i, round(p.i,2)),
                  c(n.ii, round(p.ii,2)),
                  c(n.non, round(p.non,2)),
                  c(n.cel, round(p.cel,2)),
                  c(n.control, round(p.control,2)))
colnames(pat.char) <- c('N', '%')
rownames(pat.char) <- c('Amputation', 'Death', 'Septic shock', 'Sex (male)',
                        'I', 'II', 'Non-NSTI', 'Cellulitis', 'Surgical control')
pat.char
print(paste(round(mean(data$Age),2), '?', round(sd(data$Age),2)))

################################################################################

# Predict Case type

# For unbalance correction
nrep <- 100
frac <- 0.9
n1 <- sum(data$Case_type == 'Surgical control');
n2 <- sum(data$Case_type == 'NSTI II');
n3 <- sum(data$Case_type == 'NSTI I');
n4 <- sum(data$Case_type == 'Non-NSTI');
n5 <- sum(data$Case_type == 'Celullitis');
nsize <- round(min(n1,n2,n3,n4,n5)*frac)

# Prepare the data
data$Case_type <- as.factor(data$Case_type)
attach(data)

# Initiate variables
ERROR <- NULL
important <- NULL
# Perform random forests
for(k in 1 : nrep){
  rf <- randomForest(data$Case_type ~ ., data = data[,which(colnames(data)=='C5/C5a'):ncol(data)], 
                     strata = data$Case_type, samp.size = rep(nsize, length(levels(Case_type))))
  ERROR[k] <- err.orig <- mean(rf$err.rate[,1])
  important <- cbind(important, rf$importance)
  progress(k/nrep*100)
}

# Prediction accuracy
print(100-CI(ERROR,ci = 0.95))

################################################################################

# Importance of each topological variable per type
ordered <- as.data.frame(sort(rowMeans(important), decreasing = TRUE))
ordered$var <- as.factor(rownames(ordered))
colnames(ordered)[1] <- "MeanDecreaseGini"
ordered$var <- factor(ordered$var, 
                      levels = ordered$var[order(ordered$MeanDecreaseGini)])
# Make a plot with the most important variables for prediction
number <- 20
ggplot(ordered[1:number,], aes(y = var, x = MeanDecreaseGini, color = 'Cytokine'), 
       ylab = NULL) + geom_point() + 
  scale_color_manual(values=c('orange')) + 
  theme(legend.title = element_blank()) + 
  # ggtitle(paste('Top', number, 'variables')) + 
  theme(axis.title.y=element_blank()) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0),
                                        add = c(7,1)))

################################################################################

detach(data)
