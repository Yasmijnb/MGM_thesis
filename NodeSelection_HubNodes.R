################################################################################

# Yasmijn Balder
# 18-09-2020

# Find hub nodes

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

# Hub nodes according to (Lu, et al. 2007)

hubs <- data[which(data$ClusteringCoefficient < 0.03 & data$Degree > 5),]
hubs[,c("name", "Type")]

