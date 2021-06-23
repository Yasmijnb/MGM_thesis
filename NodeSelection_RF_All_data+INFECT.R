################################################################################

# Yasmijn Balder
# 27-07-2020

# RandomForest 

################################################################################

library(randomForest)
library(Rmisc)    # For CI (confidence interval)
library(svMisc)   # For progress, not strictly necessary

################################################################################

# Load the data
data2 <- read.csv("../Provided Data/Merged/Node table - ALL+INFECT - Lambda x 2.csv")
data2$Type <- rep(0, nrow(data2))
for (sample in 1:nrow(data2)) {
  if (data2[sample,'color'] == 'orange') {
    data2[sample,'Type'] <- 'Cytokine'
  }
  if (data2[sample,'color'] == 'red') {
    data2[sample,'Type'] <- 'Human gene'
  }
  if (data2[sample,'color'] == 'green') {
    data2[sample,'Type'] <- 'Bacterial gene'
  }
  if (data2[sample,'color'] == 'cyan') {
    data2[sample,'Type'] <- 'Clinical'
  }
}
length(which(data2$color=='orange'))
length(which(data2$color=='red'))
length(which(data2$color=='green'))
length(which(data2$color=='cyan'))
data2 <- subset(data2, select = -c(frame.color,IsSingleNode,label,size,color,
                                   label.color,id,NumberOfDirectedEdges,size2,
                                   selected,SelfLoops,shape,shared.name,
                                   PartnerOfMultiEdgedNodePairs,
                                   NumberOfUndirectedEdges))

################################################################################

# For unbalance correction
nrep = 100
frac = 0.9
n1 = sum(data2$Type == 'Human gene');
n2 = sum(data2$Type == 'Bacterial gene');
nsize = round(min(n1,n2)*frac)

# Prepare the data
gene.data <- data2[which(data2$Type == 'Human gene' | data2$Type == 'Bacterial gene'),]
gene.data <- gene.data[,-which(colnames(gene.data)=='name')]
gene.data$Type <- as.factor(gene.data$Type)

# Initiate variables
ERROR = NULL
important = NULL
# Perform random forests
for(k in 1 : nrep){
  rf <- randomForest(gene.data$Type ~ ., data = gene.data, 
                     strata = gene.data$Type, samp.size = c(nsize,nsize))
  ERROR[k] = err.orig = mean(rf$err.rate[,1])
  important = cbind(important, rf$importance)
  progress(k/nrep*100)
}
# Importance of each topological variable per type
as.data.frame(rowMeans(important))
# Prediction accuracy
print(100-mean(ERROR))
print(100-CI(ERROR,ci = 0.95))

################################################################################

attach(data2)
sorted.data2 <- data2[order(Eccentricity, Radiality, AverageShortestPathLength, 
                            ClosenessCentrality, NeighborhoodConnectivity, 
                            BetweennessCentrality, TopologicalCoefficient, 
                            Stress, ClusteringCoefficient, Degree),]
detach(data2)
sorted.data2[1:20,c('name', 'Type')]
sorted.data2$name[which(sorted.data2$Type=='Clinical')][1:20]

