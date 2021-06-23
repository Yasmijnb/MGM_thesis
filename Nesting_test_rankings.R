################################################################################

# Yasmijn Balder
# 28-09-2020

# Test for nesting by looking at the order of the edge weights

################################################################################

library(arsenal)    # For comparing data frames

################################################################################

# Obtain adjacency matrices

small <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/small",
                  check.names = FALSE, row.names = 1)
big <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/big.csv",
                check.names = FALSE, row.names = 1)

# Obtain separate files (don't remove disconnected nodes)
B0 = npyLoad("B_0.npy")
Rho0 = npyLoad("Rho_0.npy")
Phi0 = npyLoad("Phi_0.npy")

# Save the small set
small <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(small) <- rownames(small) <- c(colnames(data),colnames(D))

# Save the large set
big <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(big) <- rownames(big) <- c(colnames(data),colnames(D))

################################################################################

# Rank the joint probabilities of the big data set
ranked.big <- cbind(sort(big),  # Gives the joint probability values
                    order(big)) # Gives the position of the values
ranked.big <- as.data.frame(ranked.big)
ranked.big <- ranked.big[which(ranked.big[,1]!=0),]
colnames(ranked.big) <- c('Value','Position')

ranked.big$Row <- rep(0,nrow(ranked.big))
ranked.big$Col <- rep(0,nrow(ranked.big))

for (sample in 1:nrow(ranked.big)) {
  this.value <- ranked.big$Position[sample]
  row <- ((this.value/nrow(big))%%1 * nrow(big))
  col <- ceiling(this.value/nrow(big))
  
  rowname <- rownames(big)[round(row)]
  colname <- colnames(big)[col]
  
  ranked.big$Row[sample] <- rowname
  ranked.big$Col[sample] <- colname
}

# Rank the joint probabilities of the small data set
ranked.small <- cbind(sort(small),  # Gives the joint probability values
                      order(small)) # Gives the position of the values
ranked.small <- as.data.frame(ranked.small)
ranked.small <- ranked.small[which(ranked.small[,1]!=0),]
colnames(ranked.small) <- c('Value','Position')

ranked.small$Row <- rep(0,nrow(ranked.small))
ranked.small$Col <- rep(0,nrow(ranked.small))

for (sample in 1:nrow(ranked.small)) {
  this.value <- ranked.small$Position[sample]
  row <- ((this.value/nrow(small))%%1 * nrow(small))
  col <- ceiling(this.value/nrow(small))
  
  if (row ==0) {
    row <- nrow(small)
  }
  
  rowname <- rownames(small)[round(row)]
  colname <- colnames(small)[col]
  
  ranked.small$Row[sample] <- rowname
  ranked.small$Col[sample] <- colname
}

# Keep only variables that are present in small
to_keep <- NULL
for (sample in 1:nrow(ranked.big)) {
  if (ranked.big$Row[sample] %in% c(colnames(small), "S100A8.x")) {
    to_keep <- c(to_keep, sample)
  }
  if (ranked.big$Col[sample] %in% c(colnames(small), "S100A8.x")) {
    to_keep <- c(to_keep, sample)
  }
}
ranked.big <- ranked.big[to_keep[duplicated(to_keep)],]

# Keep only even numbered rows (contains both: from A to B; from B to A)
rownames(ranked.small) <- 1:nrow(ranked.small)
ranked.small <- ranked.small[which((1:nrow(ranked.small))%% 2 == 0),]
rownames(ranked.big) <- 1:nrow(ranked.big)
ranked.big <- ranked.big[which((1:nrow(ranked.big))%% 2 == 0),]

View(ranked.small)
View(ranked.big)

################################################################################

comparison <- summary(comparedf(as.data.frame(ranked.small[,-c(1,2)]), 
                                as.data.frame(ranked.big[,-c(1,2)]), 
                                by = 'Row'))
comparison$comparison.summary.table

