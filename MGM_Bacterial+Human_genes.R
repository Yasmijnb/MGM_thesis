################################################################################

# Load all packages
library(RcppCNPy)
library(qgraph)
library(igraph)
library(RCy3)
library(missForest)

################################################################################

# Load the data
data <- read.csv("../Provided Data/Merged/Bact+Human_genes.csv",
                 check.names = FALSE)
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[,-which(colnames(data)=="RNA_Seq.data.available")]
data <- data[,-which(colnames(data)=="Most.abundant.Streptococcus.sp.")]

################################################################################

# When there is no entry, change it to NA
for (column in 1:ncol(data)) {
  # New R version:
  empty <- which(data[, column]=="")
  data[empty, column] <- NA
}

# Keep only rows which have less than a limit missing data
to_keep <- NULL
limit <- 0 # percentage of missing data that is allowed per variable
for (i in 1:nrow(data)){
  percentage <- 100*sum(is.na(data[i,]))/ncol(data)
  if (percentage<=limit) {
    to_keep <- c(to_keep, i)
  }
}
if (limit < 100) {
  data <- data[to_keep,]
}

# Keep only columns which have less than a limit missing data
# to_keep <- NULL
# limit <- 0 # percentage of missing data that is allowed per variable
# for (i in 1:ncol(data)){
#   percentage <- 100*sum(is.na(data[,i]))/nrow(data)
#   if (percentage<=limit) {
#     to_keep <- c(to_keep, i)
#   }
# }
# if (limit < 100) {
#   data <- data[,to_keep]
# }

# Remove uninteresting columns from the data
data <- subset(data, select = -c(PatientId, Day))

# Change column names for better understanding
colnames(data)[which(colnames(data) == 'Case')] <- 'City'
colnames(data)[which(colnames(data) == 'Predominant_species_humann2')] <- 'Predominant_species'

################################################################################

# Select the annotation variables and make a new data frame
Annotation <- data[,c('Location','Tissue','City',"Conclusion_micro",
                      'Classification_16S', 'Predominant_species')]

# Remove samples that have missing data points in the annotation columns
NA_samples <- NULL
for (sample in 1:nrow(Annotation)) {
  if (any(is.na(Annotation[sample,]))==TRUE) {
    NA_samples <- c(NA_samples, sample)
  }
}
if (is.null(NA_samples)==FALSE) {
  data <- data[-NA_samples,]
  Annotation <- Annotation[-NA_samples,]
}

levels <- NULL
remove_levels <- NULL
# See how many levels each variable has
for (column in 1:ncol(Annotation)) {
  level <- length(unique(Annotation[,column]))
  levels <- c(levels, level)
  if (level==1) {
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
Annotation <- Annotation[,-remove_levels]
levels <- levels[-remove_levels]

# Deselect the annotation variables from the data
data <- subset(data, select = -c(Location,Conclusion_micro,Classification_16S,
                                 Tissue,City,Predominant_species))

# Impute the missing data (default setting with 100 trees because small data-set)
# data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE)
# data_imp$OOBerror
# data <- as.data.frame(data_imp$ximp)

# # Set a limit for the amount of variables to use
# p <- 200
# # Calculate the mean stdev, and coefficient of variation of each column
# means <- colMeans(data)
# stdev <- sapply(data,sd)
# cv <- stdev/means
# # Order the genes
# data <- data[,order(cv, decreasing = TRUE)]
# # Select the top genes
# data <- data[,1:p]

# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################

# Make a discrete variable matrix
D <- NULL
# Make colnames for D
colnames_D <- NULL
# For each binary entry in the Annotation data frame
for (column in colnames(Annotation)) {
  # All of these are binary
  for(i in sort(unique(Annotation[,column]))) {
    colnames_D <- c(colnames_D, paste(column, '=', i))
    # Make all samples 0
    Dtemp <- rep(0,nrow(Annotation))
    # If the sample is 1, change it to 1
    Dtemp[Annotation[,column]==i] <- 1
    D <- rbind(D,Dtemp)
  }
}
# Transpose D so it looks the same as the original dataset 
# (samples as rows, variables as columns)
D <- as.data.frame(t(D))
# Add names to the columns
colnames(D) <- colnames_D

################################################################################

# Save matrices as npy-files in Python-script folder (here: "MGM_algorithm")
myPath <- "../Try 1/MGM_algorithm/"
setwd(myPath)
# Continuous variable matrix with rows = samples, columns = genes:
npySave("Xsc.npy", as.matrix(data))
# Discrete variable matrix with rows = samples, columns = discrete levels:
npySave("Dsc.npy", as.matrix(D))
# Save the levels vector
npySave("levelssc.npy", as.integer(levels))
# number of samples
n <- nrow(data)
# number of continuous variables
p <- ncol(data)
# number of discrete variables
q <- ncol(Annotation)
# Define lambda sequence for penalization according to (Lee and Hastie, 2015):
delta <- 1
la_seq <- delta*sqrt(log(p+q)/n)
# Save lambda sequence
npySave("lam_seqsc.npy", la_seq)
# Carry out MGM
system(paste0("python apply_functions_command_line.py",
              " --iterations=1000",           # set number of iterations to 1000
              " --oTol=1e-4",                 # set precision to 1e-4
              " --X_file=Xsc.npy",            # continuous variable matrix
              " --D_file=Dsc.npy",            # discrete variable matrix
              " --lamseq_file=lam_seqsc.npy", # lambda sequence
              " --levels_file=levelssc.npy",  # levels of discrete variables
              " --results_folder=temp"),      # output folder
       wait=TRUE)

# Move to the output folder
setwd(paste0(myPath,"temp/"))
# Matrix of edges between continuous - continuous variables:
B0 = npyLoad("B_0.npy")
# Matrix of edges between continuous - discrete variables:
Rho0 = npyLoad("Rho_0.npy")
# Matrix of edges between discrete - discrete variables:
Phi0 = npyLoad("Phi_0.npy")

################################################################################

# Combine single result matrices into one adjacency matrix:
adj_matrix <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(adj_matrix) <- rownames(adj_matrix) <- c(colnames(data),colnames(D))

# Make different groups
groups <- as.vector(c(rep("continuous",ncol(data)), rep("discrete",ncol(D))))
# Shape discrete variables nodes differently
shapes <- as.vector(c(rep("square",ncol(data)), rep("circle",ncol(D))))
# Colour discrete variables nodes differently
colours <- as.vector(c(rep("#E69F00",ncol(data)), rep("#56B4E9",ncol(D))))

# Get the column names
human_gene_data <- read.csv("../Provided Data/filt2ndHumanHgnc.csv", 
                            strip.white = TRUE, stringsAsFactors=FALSE,
                            colClasses=c('character',rep('numeric', 102)))
bact_gene_data <- read.csv("../Provided Data/Merged/Bacterial_gene_data.csv",
                           check.names = FALSE)

# Also change the shape and colour of continuous clinical parameters
for (clinical in c("Severity")) {
  colours[which(colnames(data)==clinical)] <- "#56B4E9"
}
# Colour the human and bact genes red and green
for (human in human_gene_data$Gene_name) {
  colours[which(colnames(data)==human)] <- '#CC79A7'
}
for (bact in colnames(bact_gene_data)) {
  colours[which(colnames(data)==bact)] <- '#0072B2'
}

# Remove disconnected nodes
disconnected.nodes <- which(apply(adj_matrix, 1, function(x){all(x==0)}))
if (length(disconnected.nodes)!=0) {
  adj_matrix <- adj_matrix[-disconnected.nodes,-disconnected.nodes]
  groups <- groups[-disconnected.nodes]
  shapes <- shapes[-disconnected.nodes]
  colours <- colours[-disconnected.nodes]
}

# Create a qgraph with layout options
qgraph_adj_mat <- qgraph(input=adj_matrix,
                         labels=colnames(adj_matrix),
                         groups=groups,
                         DoNotPlot=TRUE,
                         borders=FALSE,
                         palette="colorblind",
                         label.font='sans',
                         posCol="#009E73",  # colour of positive edges
                         negCol="#D55E00",  # colour of negative edges
                         color=colours,     # colour of groups
                         shape=shapes,      # shapes of groups
                         fade=FALSE,        # edge transparency based on weight
                         esize=2)

# Convert qgraph to igraph object
igraph_adj <- as.igraph(qgraph_adj_mat, attributes = TRUE)
V(igraph_adj)$name <- colnames(adj_matrix)

# Connect to cytoscape (Make sure cytoscape is opened)
cytoscapePing()
# Create the network
createNetworkFromIgraph(igraph = igraph_adj,
                        title=paste0("Full network (lambda x ", delta,')'),
                        collection="Bacterial and human genes")

################################################################################

# Make first order neighbourhood networks of the discrete variables (Cytoscape)
# How many continuous variables are left
continuous <- length(groups[groups=='continuous'])
# Turn all edges that are not connected to a discrete variable into 0
for (i in 1:continuous) {
  for (j in 1:continuous) {
    adj_matrix[i,j]<-0
  }
}

# Remove the disconnected nodes
disconnected.nodes <- which(apply(adj_matrix, 1, function(x){all(x==0)}))
if (length(disconnected.nodes)!=0) {
  adj_matrix <- adj_matrix[-disconnected.nodes,-disconnected.nodes]
  groups <- groups[-disconnected.nodes]
  shapes <- shapes[-disconnected.nodes]
  colours <- colours[-disconnected.nodes]
}

# Create a qgraph with layout options
qgraph_adj_mat <- qgraph(input=adj_matrix,
                         labels=colnames(adj_matrix),
                         groups=groups,
                         DoNotPlot=TRUE,
                         borders=FALSE,
                         label.font='sans',
                         posCol="#009E73",   # colour of positive edges
                         negCol="#D55E00",   # colour of negative edges
                         color=colours,      # colour of groups
                         shape=shapes,       # shapes of groups
                         fade=FALSE,         # edge transparency based on weight
                         esize=2)
# Convert qgraph to igraph object
igraph_adj <- as.igraph(qgraph_adj_mat, attributes = TRUE) 
V(igraph_adj)$name <- colnames(adj_matrix)

# Connect to cytoscape (Make sure cytoscape is opened)
cytoscapePing()
# Create the network
createNetworkFromIgraph(igraph = igraph_adj,
                        title=paste0('First-order neighbourhoods (lambda x ', 
                                     delta,')'), 
                        collection="Bacterial and human genes")

################################################################################

# Make first order neighbourhood networks of the discrete variables (plot in R)
# Extract the first order neighborhoods:
nei_node_graphs_adj <- make_ego_graph(graph = igraph_adj, order = 1)
# How many clinical variables are left
clinical <- which(colours=='#56B4E9')
# Plot the first order neighbourhood for each discrete variable
for (i in clinical) {
  total = names(nei_node_graphs_adj[[i]][1])
  circle = layout_in_circle(nei_node_graphs_adj[[i]], 
                            # Remove the first node from the circle
                            order=total[-which(total==colnames(adj_matrix)[i])])
  # Adjust the networks when there is a small number of nodes
  if (length(total)<4) {
    circle = layout_in_circle(nei_node_graphs_adj[[i]])
  }
  plot(nei_node_graphs_adj[[i]], main=colnames(adj_matrix)[i], layout=circle)
}

