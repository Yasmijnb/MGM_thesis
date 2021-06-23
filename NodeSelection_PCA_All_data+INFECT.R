# PCA on node table

library(ggplot2)    # For function 'autoplot' to plot the PCA
# library(devtools)
# install_github('sinhrks/ggfortify')
library(ggfortify)  # For autoplot options?

################################################################################

# Load the data
data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Node table - ALL+INFECT - Delta = 2.csv")
data$Type <- rep(0, nrow(data))
for (sample in 1:nrow(data)) {
  if (data[sample,'color'] == 'orange') {
    data[sample,'Type'] <- 'Cytokine'
  }
  if (data[sample,'color'] == 'red') {
    data[sample,'Type'] <- 'Human gene'
  }
  if (data[sample,'color'] == 'green') {
    data[sample,'Type'] <- 'Bacterial gene'
  }
  if (data[sample,'color'] == 'cyan') {
    data[sample,'Type'] <- 'Clinical'
  }
}
length(which(data$color=='orange'))
length(which(data$color=='red'))
length(which(data$color=='green'))
length(which(data$color=='cyan'))
data <- subset(data, select = -c(frame.color,IsSingleNode,label,size,color,
                                   label.color,id,NumberOfDirectedEdges,size2,
                                   selected,SelfLoops,shape,shared.name,
                                   PartnerOfMultiEdgedNodePairs,
                                   NumberOfUndirectedEdges))

################################################################################

# Make a PCA for lambda x 2
pca <- prcomp(data[,-c(7,12)], center = TRUE, scale = TRUE)

summary(pca)$importance[2:3,]

# Plot the first and second PCs
autoplot(pca, data = data, colour = 'Type', loadings = TRUE, legend = TRUE, 
         x = 1, y = 2, loadings.label = TRUE) + 
  # labs(title="Lambda x 2") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title = element_blank()) + 
  scale_colour_manual(values = c('dark green','cyan','orange','red'))

# Proportion of Variance Explained
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
plot(pve, type = "o", ylab = "PVE", xlab = "Principal Component", col = "blue")
plot(cumsum(pve), type = "o", ylab = "Cumulative PVE", 
     xlab = "Principal Component", col = "brown3")

# See which loadings attribute to the first PC
PC1.loadings <- pca$rotation[,1]
# To get rid of negative and positive values, square the numbers
PC1.loadings <- (PC1.loadings)^2
# Sort the samples by most attribution to PC1 
PC1.loadings <- sort(PC1.loadings, decreasing = TRUE)
# View the top attributing loadings
PC1.loadings

# # See which loadings attribute to the second PC
# PC2.loadings <- pca$rotation[,2]
# # To get rid of negative and positive values, square the numbers
# PC2.loadings <- (PC2.loadings)^2
# # Sort the samples by most attribution to PC1 
# PC2.loadings <- sort(PC2.loadings, decreasing = TRUE)
# # View the top attributing loadings
# PC2.loadings

# View(PC1.loadings)
# View(PC2.loadings)

################################################################################

# Sort the nodes by the highest loadings
attach(data)
sorted.data <- data[order(-AverageShortestPathLength, Radiality, 
                            -Eccentricity, ClosenessCentrality, 
                            -TopologicalCoefficient, -ClusteringCoefficient, 
                            -NeighborhoodConnectivity, -Degree, 
                            BetweennessCentrality, Stress),]
detach(data)
sorted.data$name[1:20]
sorted.data$name[which(sorted.data$Type=='Clinical')][1:10]