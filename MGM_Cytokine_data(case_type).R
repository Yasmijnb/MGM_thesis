################################################################################

# MGM with cytokines, including controls

################################################################################

# Load all packages
library(RcppCNPy)
library(qgraph)
library(igraph)
library(RCy3)

################################################################################

# Load the data
library("readxl")
data <- read_excel("../Provided Data/Cytokines_3Dec2019.xlsx")
data <- as.data.frame(data)
data <- data[,-which(colnames(data)=="ID_Day")]
data <- data[,-which(colnames(data)=="ID_day")]
data <- data[data$day == 0,]
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
for (i in 1:ncol(data)){
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
                 'SOFA6','SOFA7','SOFA1','SAPSII', 'DeathOnly')) {
  data <- data[,-which(colnames(data)==column)]
}

################################################################################

# Remove some rows
data <- data[-which(data$BMI=='NA'),]

# For nesting test
# data <- data[data$Case_type == 'NSTI I' | data$Case_type == 'NSTI II',]

################################################################################

# Separate continuous and discrete variables
Annotation <- data[,c("Case_type","Sex","Comorbidity","Septicshock","IVIG",
                      "Clindamycin")]
data <- subset(data, select = -c(Case_type,Sex,Comorbidity,Septicshock,IVIG,
                                 Clindamycin))

levels <- NULL
remove_levels <- NULL
# See how many levels each variable has
for (column in 1:ncol(Annotation)) {
  level <- length(unique(Annotation[,column]))
  levels <- c(levels, level)
  if (level==1) {
    # print(colnames(Annotation)[column])
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
Annotation <- Annotation[,-remove_levels]
levels <- levels[-remove_levels]

class(data$BMI) <- 'numeric'

# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################

# Make a discrete variable matrix
D <- NULL
# Make colnames for D
colnames_D <- NULL
# For each binary entry in the Annotation data frame
# Make a discrete variable matrix
D <- NULL
colnames_D <- NULL
# For each discrete entry in the Annotation data frame
for (colnum in 1:ncol(Annotation)) {
  # If the variable has names as levels
  # if (is.null(levels(Annotation[,colnum]))==FALSE) {
  if (class(Annotation[,colnum])=='character') {
    # print(paste(colnum, 'character', levels[colnum]))
    for (i in sort(unique(Annotation[,colnum]))) {
      # Make all samples 0
      Dtemp <- rep(0,nrow(Annotation))
      # If the sample is this level, change it to 1
      Dtemp[Annotation[,colnum]==i] <- 1
      D <- rbind(D,Dtemp)
      # Add the colname
      colnames_D <- c(colnames_D, paste(colnames(Annotation)[colnum],'=',i))
  }} else {
    # If the variable has numbers as levels
    for (i in sort(unique(Annotation[,colnum]))) {
      # Start at the minimum level
      # mini <- min(Annotation[,colnum])
      # i <- i+mini
      # Make all samples 0
      Dtemp <- rep(0,nrow(Annotation))
      # If the sample is this level, change it to 1
      Dtemp[as.integer(Annotation[,colnum])==i] <- 1
      D <- rbind(D,Dtemp)
      # Add the colname
      colnames_D <- c(colnames_D, paste(colnames(Annotation)[colnum],'=',i))
      }
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
delta <- 2
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
shapes <- as.vector(c(rep("circle",ncol(data)), rep("square",ncol(D))))
# Colour discrete variables nodes differently
colours <- as.vector(c(rep("#E69F00",ncol(data)), rep("#56B4E9",ncol(D))))

# Also change the shape and colour of continuous clinical parameters
for (clinical in c("Age","BMI")) {
  colours[which(colnames(data)==clinical)] <- "#56B4E9"
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
                        collection="Cytokines")

################################################################################

# Make first order neighbourhood networks of the discrete variables (Cytoscape)
# How many continuous variables are left
continuous <- which(colours!='#56B4E9')
# Turn all edges that are not connected to a clinical variable into 0
for (i in continuous) {
  for (j in continuous) {
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
                        title=paste0('First-order neighbourhoods of clinical parameters (lambda x ',
                                     delta,') (limit = ', limit, ')'),
                        collection="Cytokines")

# # Extract the first order neighborhoods:
nei_node_graphs_adj <- make_ego_graph(graph = igraph_adj, order = 1)
# How many clinical variables are left
clinical <- which(colours=='cyan')
# Plot the first order neighbourhood for each discrete variable
for (i in clinical) {
  total <- names(nei_node_graphs_adj[[i]][1])
  circle <- layout_in_circle(nei_node_graphs_adj[[i]],
                             # Remove the first node from the circle
                             order=total[-which(total==colnames(adj_matrix)[i])])
  # Adjust the networks when there is a small number of nodes
  if (length(total)<4) {
    # Make sure that the vertices are on one line
    for (coord in 1:nrow(circle)) {
      if (circle[coord,1]==-1) {
        circle[coord,2] <- 0
      }
    }
    # Reorder so that a vertical instead of a horizontal network is created
    circle <- cbind(circle[,2], circle[,1])
  }
  # Plot the network
  plot(nei_node_graphs_adj[[i]], main=colnames(adj_matrix)[i], layout=circle)
}



