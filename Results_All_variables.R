
# Make table as in appendix of helena's paper
# Run the preprocessing, but the MGM is not necessary. Don't scale the data.

################################################################################

# Initiate columns
Variablename <- NULL
Type <- NULL
Category <- NULL
Continuous <- NULL
Discrete <- NULL

################################################################################

# Prepare

colours <- as.vector(rep("Clinical",ncol(data)))
bact_gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Bacterial_gene_data.csv",
                           check.names = FALSE)
human_gene_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/filt2ndHumanHgnc.csv", 
                            strip.white = TRUE, stringsAsFactors=FALSE,
                            colClasses=c('character',rep('numeric', 102)))
cytokine_data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Final_Data_Imputed.csv")
for (human in human_gene_data$Gene_name) {
  colours[which(colnames(data)==human)] <- 'human gene'
}
for (bact in colnames(bact_gene_data)) {
  colours[which(colnames(data)==bact)] <- 'bacterial gene'
}
for (cyto in colnames(cytokine_data)) {
  colours[which(colnames(data)==cyto)] <- 'cytokine'
}
for (clinical in c("Severity", "SAPSII", 'Age', 'SOFA1', 'BMI', 'IVIG', 'Comorbidity')) {
  colours[which(colnames(data)==clinical)] <- "clinical"
}
for (odd in c('S100A9', 'ICAM1', 'S100A8')) {
  loc.x <- which(colnames(data)==paste0(odd,'.x'))
  loc.y <- which(colnames(data)==paste0(odd,'.y'))
  colours[loc.x] <- 'cytokine'
  colours[loc.y] <- 'human gene'
}

########################################

# Continuous
for (colnum in 1:ncol(data)) {
  Variablename <- c(Variablename, colnames(data)[colnum])
  Type <- c(Type, 'continous')
  Category <- c(Category, colours[colnum])
  m <- mean(data[,colnum])
  s <- sd(data[,colnum])
  Continuous <- c(Continuous, as.character(paste(signif(m,6), signif(s,6), sep = ' ± ')))
  Discrete <- c(Discrete, '-')
}

################################################################################

# Discrete
for (colnum in 1:ncol(Annotation)) {
  name <- str_remove_all(string = colnames(Annotation)[colnum], pattern = '\\((\\w+\\=[\\w\\s\\-\\.]+[\\,\\s\\)]+)+')
  name <- trimws(name, which = 'both')
  Variablename <- c(Variablename, name)
  Type <- c(Type, 'discrete')
  Category <- c(Category, 'clinical')
  Continuous <- c(Continuous, '-')
  uni <- sort(unique(Annotation[,colnum]))
  string <- NULL
  for (i in uni) {
    string <- paste0(string, i, sep = '; ')
  }
  Discrete <- c(Discrete, string)
}

################################################################################

# Combine into table
overview <- cbind(Variablename, Type, Category, Continuous, Discrete)
overview <- as.data.frame(overview)

# View(overview)

# Save table
myPath <- "~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/"
setwd(myPath)
write.csv(overview, 'Overview of variables.csv')
