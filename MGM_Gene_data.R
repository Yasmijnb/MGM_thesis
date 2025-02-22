################################################################################

# Load all packages
library(RcppCNPy)
library(qgraph)
library(igraph)
library(RCy3)
# library(impute)
library(missForest)

################################################################################

# Load the data
data <- read.csv("../Provided Data/common_Gene.txt"
                 , sep = '\t', skip = 1)
data <- data[,-which(colnames(data)=='ID')]

# Select the annotation variables and make a new data frame
Annotation <- data[,c('Type_num','Sex_y','Comorbidity_y',"Amputation_y",'Death',
                      "Septicshock_at_Baseline","IVIG_y","Clindamycin_y")]

colnames(Annotation)<-c('NSTI type','Sex','Comorbidity',"Amputation",'Death',
                        "Septicshock at baseline","IVIG","Clindamycin")

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
data <- subset(data, select = -c(ID,Type_num,Sex_y,Comorbidity_y,Amputation_y,
                                 Septicshock_at_Baseline,Clindamycin_y,IVIG_y,
                                 Death))

# See how many levels each variable has
for (column in 1:ncol(data)) {
  level <- length(unique(data[,column]))
  if (level==1) {
    remove_levels <- c(remove_levels,column)
  }
}
# Remove columns with only one level
data <- data[,-remove_levels]

################################################################################

# Print some patient characteristics
n.male <- length(which(Annotation$Sex == 1))
p.male <- n.male*100/nrow(data)
n.ss <- length(which(Annotation$`Septicshock at baseline` == 1))
p.ss <- n.ss*100/nrow(data)
n.amp <- length(which(Annotation$Amputation == 1))
p.amp <- n.amp*100/nrow(data)
n.death <- length(which(Annotation$Death == 1))
p.death <- n.death*100/nrow(data)
pat.char <- rbind(c(n.male, p.male), 
                  c(n.ss, p.ss),
                  c(n.amp, p.amp),
                  c(n.death, p.death))
colnames(pat.char) <- c('Median / N', 'Quantile / %')
rownames(pat.char) <- c('Sex (male)', 'Septic shock', 'Amputation', 'Death')
pat.char
print(paste(round(mean(data$Age),2), '?', round(sd(data$Age),2)))

################################################################################

# Impute the missing data (default setting with 100 trees because small dataset)
data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE)
data_imp$OOBerror
data <- as.data.frame(data_imp$ximp)

# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################

# Make a discrete variable matrix
D <- NULL
# For each binary entry in the Annotation data frame
for (column in colnames(Annotation)) {
  # All of these are binary
  for(i in c(0:1)){
    # Make all samples 0
    Dtemp <- rep(0,nrow(Annotation))
    # If the sample is 1, change it to 1
    Dtemp[as.integer(Annotation[,column])==i] <- 1
    D <- rbind(D,Dtemp)
  }
}
# Transpose D so it looks the same as the original dataset 
# (samples as rows, variables as columns)
D <- as.data.frame(t(D))
# Make colnames for D
colnames_D <- NULL
# For each discrete variable
for (i in colnames(Annotation)) {
  # Make two colnames, =0 and =1 (without spaces)
  colnames_D <-c(colnames_D, as.vector('baseline'), as.vector(paste0(i)))
}
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
shapes <- as.vector(c(rep("circle",ncol(data)), rep("square",ncol(D))))
# Colour discrete variables nodes differently
colours <- as.vector(c(rep("#CC79A7",ncol(data)), rep("#56B4E9",ncol(D))))

# Also change the shape and colour of continuous clinical parameters
for (clinical in c("Age_y","BMI_y","SAPSII","SOFA1")) {
  colours[which(colnames(data)==clinical)] <- "#56B4E9"
}

# Remove diconnected nodes
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

# Connect to cytoscape (Make sure cytoscape is opened)
cytoscapePing()
# Create the network
createNetworkFromIgraph(igraph = igraph_adj,
                        title=paste0("Full network (lambda x ", delta,')'),
                        collection="common_Gene")

################################################################################

# Make first order neighbourhood networks of the discrete variables (Cytoscape)
# How many cytokine variables are left
cytokines <- which(colours=='#CC79A7')
# Turn all edges that are not connected to a clinical variable into 0
for (i in cytokines) {
  for (j in cytokines) {
    adj_matrix[i,j]<-0
  }
}

# Remove the diconnected nodes
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
                        collection="common_Gene")

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
    # Make sure that the vertices are on one line
    for (coord in 1:nrow(circle)) {
      if (circle[coord,1]==-1) {
        circle[coord,2] <- 0
      }
    }
    # Reorder so that a vertical instead of a horizontal network is created
    circle <- cbind(circle[,2], circle[,1])
  }
  plot(nei_node_graphs_adj[[i]], main=colnames(adj_matrix)[i], layout=circle)
}
