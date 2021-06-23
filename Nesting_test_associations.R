################################################################################

# Yasmijn Balder
# 14-09-2020

# Test for nesting by looking at the joint probability conclusions

################################################################################
library(arsenal)    # For comparing data frames

################################################################################

# Obtain adjacency matrices

small <- read.csv("../Provided Data/Merged/small", check.names = FALSE, row.names = 1)
big <- read.csv("../Provided Data/Merged/big.csv", check.names = FALSE, row.names = 1)

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
# 
# # Compare adjacency matrices on edge weights
# 
# comparedf(as.data.frame(small), as.data.frame(big), by = 'row.names')
# comparison <- summary(comparedf(as.data.frame(small), as.data.frame(big), 
#                                 by = 'row.names'))
# comparison$comparison.summary.table
# 
################################################################################

# Convert any positive number to 1, and any negative number to -1

# For the small data set
for (colnum in 1:ncol(small)) {
  if (colnum == 1) {
    b.small <- as.data.frame(rep(0, nrow(small)))
    b.small[small[,colnum] > 0,1] <- 1
    b.small[small[,colnum] < 0,1] <- -1
  } 
  else{
    b.small$new <- 0
    b.small$new[small[,colnum] > 0] <- 1
    b.small$new[small[,colnum] < 0] <- -1 
  }
  colnames(b.small)[ncol(b.small)] <- colnames(small)[colnum]
}
rownames(b.small) <- colnames(b.small)

# For the large data set
for (colnum in 1:ncol(big)) {
  if (colnum == 1) {
    b.big <- as.data.frame(rep(0, nrow(big)))
    b.big[big[,colnum] > 0,1] <- 1
    b.big[big[,colnum] < 0,1] <- -1
  } 
  else{
    b.big$new <- 0
    b.big$new[big[,colnum] > 0] <- 1
    b.big$new[big[,colnum] < 0] <- -1 
  }
  colnames(b.big)[ncol(b.big)] <- colnames(big)[colnum]
}
rownames(b.big) <- colnames(b.big)

################################################################################

# Compare adjacency matrices on conclusion
comparison <- summary(comparedf(as.data.frame(b.small), as.data.frame(b.big), 
              by = 'row.names'))
comparison$comparison.summary.table
