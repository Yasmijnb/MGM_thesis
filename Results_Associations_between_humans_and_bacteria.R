################################################################################

# Yasmijn Balder
# 10-11-2020

# Associations between human and bacteria

################################################################################

# Load all packages
library("readxl")           # For reading the cytokine excel file
library(ggplot2)            # For making nice plots

################################################################################

# Run the MGM of the preferred data set first

# Obtain all edges (don't remove disconnected nodes)
myPath <- "../Try 1/MGM_algorithm/"
setwd(paste0(myPath,"temp/"))
B0 = npyLoad("B_0.npy")
Rho0 = npyLoad("Rho_0.npy")
Phi0 = npyLoad("Phi_0.npy")

# Save the adjacency matrix
adj_matrix <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(adj_matrix) <- rownames(adj_matrix) <- c(colnames(data),colnames(D))

################################################################################

# Rank the joint probabilities of the adj_matrix data set
ranked.adj_matrix <- cbind(sort(adj_matrix, decreasing = TRUE),  # Gives the joint probability values
                           order(adj_matrix, decreasing = TRUE)) # Gives the position of the values
ranked.adj_matrix <- as.data.frame(ranked.adj_matrix)
ranked.adj_matrix <- ranked.adj_matrix[which(ranked.adj_matrix[,1]!=0),]
colnames(ranked.adj_matrix) <- c('Value','Position')

ranked.adj_matrix$Row <- rep(0,nrow(ranked.adj_matrix))
ranked.adj_matrix$Col <- rep(0,nrow(ranked.adj_matrix))

for (sample in 1:nrow(ranked.adj_matrix)) {
  this.value <- ranked.adj_matrix$Position[sample]
  row <- ((this.value/nrow(adj_matrix))%%1 * nrow(adj_matrix))
  col <- ceiling(this.value/nrow(adj_matrix))
  
  if (col ==0) {
    col <- nrow(adj_matrix)
  }
  
  if (row ==0) {
    row <- nrow(adj_matrix)
  }
  
  rowname <- rownames(adj_matrix)[round(row)]
  colname <- colnames(adj_matrix)[col]
  
  ranked.adj_matrix$Row[sample] <- rowname
  ranked.adj_matrix$Col[sample] <- colname
}

# Keep only even numbered rows (contains both: from A to B; from B to A)
rownames(ranked.adj_matrix) <- 1:nrow(ranked.adj_matrix)
# Make column with association
ranked.adj_matrix$Association <- rep(NA, nrow(ranked.adj_matrix))
for (sample in 1:nrow(ranked.adj_matrix)) {
  row <- ranked.adj_matrix$Row[sample]
  col <- ranked.adj_matrix$Col[sample]
  if (sample%% 2 == 0) {
    ranked.adj_matrix$Association[sample] <- paste0(row,'\n',col)
  } else {
    ranked.adj_matrix$Association[sample] <- paste0(col,'\n',row)
  }
}

################################################################################

# Get the names
bact_gene_data <- read.csv("../Provided Data/Merged/Bacterial_gene_data.csv", 
                           check.names = FALSE)
human_gene_data <- read.csv("../Provided Data/filt2ndHumanHgnc.csv", 
                            strip.white = TRUE, stringsAsFactors=FALSE,
                            colClasses=c('character',rep('numeric', 102)))
cytokine_data <- read_excel("../Provided Data/Cytokines_3Dec2019.xlsx")

# Make a selection 
to_keep <- NULL
all.human <- c(human_gene_data$Gene_name,'S100A9.y','ICAM1.y','S100A8.y')
all.cytos <- c(colnames(cytokine_data)[28:66],'S100A9.x','ICAM1.x','S100A8.x')
for (sample in 1:nrow(ranked.adj_matrix)) {
  # Bacterial gene with human gene
  if (ranked.adj_matrix$Row[sample] %in% colnames(bact_gene_data)[1:1580] &&
      ranked.adj_matrix$Col[sample] %in% all.human) {
    to_keep <- c(to_keep, sample)
  }
  # Human gene with bacterial gene
  if (ranked.adj_matrix$Col[sample] %in% colnames(bact_gene_data)[1:1580] &&
      ranked.adj_matrix$Row[sample] %in% all.human) {
    to_keep <- c(to_keep, sample)
  }
  # Bacterial gene with cytokine
  if (ranked.adj_matrix$Row[sample] %in% colnames(bact_gene_data)[1:1580] &&
      ranked.adj_matrix$Col[sample] %in% all.cytos) {
    to_keep <- c(to_keep, sample)
  }
  # Cytokine with bacterial gene
  if (ranked.adj_matrix$Col[sample] %in% colnames(bact_gene_data)[1:1580] &&
      ranked.adj_matrix$Row[sample] %in% all.cytos) {
    to_keep <- c(to_keep, sample)
  }
}
ranked.adj_matrix <- ranked.adj_matrix[to_keep,]

################################################################################

# Add the rank
ranked.adj_matrix$Rank <- 1:nrow(ranked.adj_matrix)
all.assoc <- unique(ranked.adj_matrix$Association)
for (sample in 1:nrow(ranked.adj_matrix)) {
  this.assoc <- ranked.adj_matrix$Association[sample]
  ranked.adj_matrix$Rank[sample] <- which(all.assoc==this.assoc)
}

# Add the variable type
ranked.adj_matrix$Type <- rep("Clinical", nrow(ranked.adj_matrix))
for (bact in colnames(bact_gene_data)) {
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==bact)] <- 'Bacterial gene'
}
for (human in human_gene_data$Gene_name) {
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==human)] <- 'Human gene'
}
for (cyto in colnames(cytokine_data)) {
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==cyto)] <- 'Cytokine'
}
for (clinical in c("Severity", "SAPSII", 'Age', 'SOFA1', 'BMI', 'IVIG', 'Comorbidity')) {
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==clinical)] <- "Clinical"
}
for (odd in c('S100A9', 'ICAM1', 'S100A8')) {
  cyto <- paste0(odd,'.x')
  human <- paste0(odd,'.y')
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==cyto)] <- 'Cytokine'
  ranked.adj_matrix$Type[which(ranked.adj_matrix$Row==human)] <- 'Human gene'
}

################################################################################

View(ranked.adj_matrix)

# Make plot
plot.colours <- c("#0072B2","#E69F00","#CC79A7")

ggplot(ranked.adj_matrix, aes(fill = Type, 
                              x = reorder(Association, -Rank), y = Rank)) +
  # Make a grouped bar plot
  geom_col(position="dodge", width = 0.5) +
  # Use these colours for the bars
  scale_fill_manual(values=plot.colours) +
  # Use horizontal bars
  coord_flip() +
  # Change the y-axis name
  labs(x = "Association") +
  # Don't include the legend title
  theme(legend.title = element_blank()) + 
  # Make a white background
  theme_bw()

################################################################################
