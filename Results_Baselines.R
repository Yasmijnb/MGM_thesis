################################################################################

# Yasmijn Balder
# 24-06-2020

# Makes a dataframe with the baselines for each discrete variable
# Run MGM code first

################################################################################

# Create an empty dataframe and give it column names
baselines <- as.data.frame(rep(0, ncol(Annotation)))
colnames(baselines) <- 'Variable'
baselines$Baseline <- rep(0, ncol(Annotation))

# Fill the dataframe with the variables and corresponding baselines
for (i in 1:ncol(Annotation)) {
  name <- str_remove_all(string = colnames(Annotation)[i], pattern = '\\((\\w+\\=[\\w\\s\\-\\.]+[\\,\\s\\)]+)+')
  name <- trimws(name, which = 'both')
  it <- sum(levels[1:i])
  lev <- colnames(D)[it - levels[i] + 1]
  splitlev <- strsplit(lev, split = '[=]')[[1]]
  for (split in 1:length(splitlev)) {
    splitlev[split] <- trimws(splitlev[split])
  }
  lev <- splitlev[2]
  baselines[i,] <- c(name, lev)
}

View(baselines)

# Save the table as a csv file
myPath <- "~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/"
setwd(myPath)
write.csv(baselines, 'Baselines_All+INFECT.csv')

# For clinical data
# for (i in c(1,2,78,114,115,122,126,127,131,134:146,148,149)) {
#   print(baselines$Baseline[i])
# }