###############################################################################

# Load all packages
library(RcppCNPy)           # Used to save and load numpy files
library(qgraph)             # Used to create the network
library(igraph)             # Used to create the network
library(RCy3)               # Used to open the network in cytoscape
library(missForest)         # Used to impute missing data
library(stringr)            # Used to edit the names
library("readxl")           # For reading the cytokine excel file

################################################################################

# Load the data
data <- read.csv("../Provided Data/Merged/All_data+INFECT.csv",
                      check.names = FALSE, row.names = 1)
data <- data[,-which(colnames(data)=='Patient ID')]

to_remove <- NULL
for (sample in 1:nrow(data)) {
  if (is.na(data$`C5/C5a`[sample])==TRUE) {
    to_remove <- c(to_remove, sample)
  }
  if (is.na(data$CFLAR[sample])==TRUE) {
    to_remove <- c(to_remove, sample)
  }
}
data <- data[-unique(to_remove),]

################################################################################

# Keep only columns which have less than a limit missing data
to_keep <- NULL
limit <- 0 # percentage of missing data that is allowed per variable
for (i in 1:ncol(data)){
  percentage <- 100*sum(is.na(data[,i]))/nrow(data)
  if (percentage<=limit) {
    # print(percentage)
    to_keep <- c(to_keep, i)
  }
}
if (limit < 100) {
  data <- data[,to_keep]
}

################################################################################

# Assign a column to annotation based on the name
discrete_columns <- grepl('=', colnames(data))
# Add these columns to the Annotation matrix
Annotation <- data[,which(discrete_columns==TRUE)]
# Remove these columns from the data matrix
data <- data[,-which(discrete_columns==TRUE)] 
# Select factors to the Annotation 
for (column in colnames(data)) {
  if (class(data[,column])=="character") {
    if (is.na(match(column, colnames(data)))==FALSE) {
      # Add this column to the Annotation matrix
      Annotation <- cbind(Annotation, data[,column])
      colnames(Annotation)[ncol(Annotation)] <- column
      # Remove this column from the data matrix
      data <- data[,-which(colnames(data)==column)]
    }
  }
}
# Select left over columns
for (column in c('Alcohol Consumption','Mechanical ventilation Baseline',
                 "Need For reconstructive Surgery ((if not in hospital day 90)",
                 'NSTI type', 'Sex', 'Comorbidity', "amputation", "Septicshock", 
                 "IVIG", "Clindamycin", 'Death',"Other Antibiotics Day 1",
                 "Other Antibiotics Given Before ICU Admission", 
                 "Other Antibiotics Day 2","Other Antibiotics Day 3", 
                 "Other Antibiotics Day 4", "Other Antibiotics Day 5", 
                 "Other Antibiotics Day 6", "Other Antibiotics Day 7", 
                 "Bodypart affected at arrival - conclusion", 
                 "Microbiological species found in tissue sample", 
                 "Microbiological agens - conclusion","LRINEC Risk of NSTI",
                 "LRINEC risk of NSTI dichotomous",
                 "Amputation of Limb (Time: from diagnosis to ICU day 7)")) {
  if (is.na(match(column, colnames(data)))==FALSE) {
    # Add this column to the Annotation matrix
    Annotation <- cbind(Annotation, data[,column])
    colnames(Annotation)[ncol(Annotation)] <- column
    # Remove this column from the data matrix
    data <- data[,-which(colnames(data)==column)]
  }
}

# Remove samples that have missing data points in the annotation columns
NA_samples <- NULL
for (sample in 1:nrow(Annotation)) {
  if (any(is.na(Annotation[sample,]))==TRUE) {
    # print(rownames(Annotation)[sample])
    NA_samples <- c(NA_samples, sample)
  }
} # When no limit is set, this step removes all samples
if (is.null(NA_samples)==FALSE) {
  data <- data[-NA_samples,]
  Annotation <- Annotation[-NA_samples,]
}

# See how many level each discrete variable has
levels <- NULL
remove_level <- NULL
for (column in 1:ncol(Annotation)) {
  level <- length(unique(Annotation[,column]))
  levels <- c(levels, level)
  if (level==1) {
    remove_level <- c(remove_level, column)
  }
}
if (is.null(remove_level)==FALSE) {
  Annotation <- Annotation[,-remove_level]
  levels <- levels[-remove_level]
}

# Remove continuous variables with only one level
remove_level <- NULL
for (column in 1:ncol(data)) {
  level <- length(unique(data[,column]))
  if (level==1) {
    remove_level <- c(remove_level, column)
  }
}
if (is.null(remove_level)==FALSE) {
  data <- data[,-remove_level]
}

################################################################################

# Print some patient characteristics
n.male <- length(which(Annotation$`Sex (0=female, 1=male)` == 1))
p.male <- n.male*100/nrow(data)
n.ss <- length(which(Annotation$`Septic shock Baseline (1=yes, 0=no)` == 1))
p.ss <- n.ss*100/nrow(data)
n.amp <- length(which(Annotation$amputation == 1))
p.amp <- n.amp*100/nrow(data)
n.death <- length(which(Annotation$`Death (0=no, 1=yes)` == 1))
p.death <- n.death*100/nrow(data)
n.mono <- length(which(Annotation$Case_type=='NSTI II'))
p.mono <- n.mono*100/nrow(Annotation)
n.poly <- length(which(Annotation$Case_type=='NSTI I'))
p.poly <- n.poly*100/nrow(Annotation)
pat.char <- rbind(c(n.amp, round(p.amp,2)),
                  c(n.death, round(p.death,2)),
                  c(n.ss, round(p.ss,2)),
                  c(n.male, round(p.male,2)),
                  c(n.mono, round(p.mono,2)),
                  c(n.poly, round(p.poly,2)))
colnames(pat.char) <- c('N', '%')
rownames(pat.char) <- c('Amputation', 'Death', 'Septic shock', 'Sex (male)',
                        'Mono', 'Poly')
pat.char
print(paste(round(mean(data$`Age (years)`),2), '?', round(sd(data$`Age (years)`),2)))

################################################################################

# Impute the missing data
if (limit > 0) {
  data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE, ntree = 100)
  data_imp$OOBerror
  data <- as.data.frame(data_imp$ximp)
}

# Remove continuous variables with only one level
remove_level <- NULL
for (column in 1:ncol(data)) {
  level <- length(unique(data[,column]))
  if (level==1) {
    remove_level <- c(remove_level, column)
  }
}
# If the variable has only 1 level, remove it from the data
if (is.null(remove_level)==FALSE) {
  data <- data[,-remove_level]
}

# Scale the continuous variables
data <- as.data.frame(scale(data))

################################################################################
## MAKE THE D MATRIX FROM THE ANNOTATION MATRIX

# Make a discrete variable matrix
D <- NULL
colnames_D <- NULL
# For each discrete entry in the Annotation data frame
for (colnum in 1:ncol(Annotation)) {
  name <- str_remove_all(string = colnames(Annotation)[colnum], pattern = '\\((\\w+\\=[\\w\\s\\-\\.]+[\\,\\s\\)]+)+')
  name <- trimws(name, which = 'both')
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
      # if (levels[colnum]==2) {
      #   colnames_D <- c(colnames_D,name)
      # } else 
      if (grepl('=',colnames(Annotation)[colnum])==TRUE) {
        splitted<-strsplit(colnames(Annotation)[colnum], split = '[=]|[(]|[)]|[,]')[[1]]
        for (split in 1:length(splitted)) {
          splitted[split] <- trimws(splitted[split])
        }
        level<-splitted[which(splitted==i)+1]
        colnames_D <- c(colnames_D,paste(name,'=',level)) 
      } else { # if (levels[colnum]!=2) {
        colnames_D <- c(colnames_D,paste(colnames(Annotation)[colnum],'=',i))}
      # print(colnames_D[length(colnames_D)])
    }
    
  } else {
    # print(paste(colnum,'integer', levels[colnum]))
    # If the variable has numbers as levels
    # for(i in 0:(levels[colnum]-1)){
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
      # if (levels[colnum]==2) {
      #   colnames_D <- c(colnames_D,name)
      # } else 
      if (grepl('=',colnames(Annotation)[colnum])==TRUE) {
        splitted<-strsplit(colnames(Annotation)[colnum], split = '[=]|[(]|[)]|[,]')[[1]]
        for (split in 1:length(splitted)) {
          splitted[split] <- trimws(splitted[split])
        }
        # if ('no' %in% splitted) {
        #   splitted <- c(splitted, 2, 'unknown')
        # }
        level<-splitted[which(splitted==i)+1]
        if (level == 'yes'){#} || level == 'yes') {
          colnames_D <- c(colnames_D,name)
        }
        if (level != 'yes'){#&& level != 'yes') {
          colnames_D <- c(colnames_D,paste(name,'=',level))
        }
      } else { # if (levels[colnum]!=2)
        if (i == 'yes') {
          colnames_D <- c(colnames_D,name)
          }
      else {
          colnames_D <- c(colnames_D,paste(colnames(Annotation)[colnum],'=',i))}
      }
      # print(colnames_D[length(colnames_D)])
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
colours <- as.vector(c(rep("#56B4E9",ncol(data)), rep("#56B4E9",ncol(D))))

# Get the names
bact_gene_data <- read.csv("../Provided Data/Merged/Bact+Human_genes.csv",
                           check.names = FALSE)
human_gene_data <- read.csv("../Provided Data/filt2ndHumanHgnc.csv", 
                            strip.white = TRUE, stringsAsFactors=FALSE,
                            colClasses=c('character',rep('numeric', 102)))
cytokine_data <- read_excel("../Provided Data/Cytokines_3Dec2019.xlsx")

# Colour the human and bact genes
for (bact in colnames(bact_gene_data)) {
  colours[which(colnames(adj_matrix)==bact)] <- '#0072B2'
}
for (human in human_gene_data$Gene_name) {
  colours[which(colnames(adj_matrix)==human)] <- '#CC79A7'
}
for (cyto in colnames(cytokine_data)) {
  colours[which(colnames(adj_matrix)==cyto)] <- '#E69F00'
}
# Change the shape and colour of continuous clinical parameters
for (clinical in c("Severity", "SAPSII", 'Age', 'SOFA1', 'BMI', 'IVIG', 'Comorbidity')) {
  colours[which(colnames(adj_matrix)==clinical)] <- "#56B4E9"
}
# Some names are in both the cytokines and the genes
for (odd in c('S100A9', 'ICAM1', 'S100A8')) {
  loc.x <- which(colnames(data)==paste0(odd,'.x'))
  loc.y <- which(colnames(data)==paste0(odd,'.y'))
  colours[loc.x] <- '#E69F00'
  colours[loc.y] <- '#CC79A7'
}

n
length(which(colours=='#CC79A7'))
length(which(colours=='#0072B2'))
length(which(colours=='#E69F00'))
p + q - length(which(colours=='#CC79A7')) - length(which(colours=='#0072B2')) - length(which(colours=='#E69F00'))
q

# Remove disconnected nodes
disconnected.nodes <- which(apply(adj_matrix, 1, function(x){all(x==0)}))
if (length(disconnected.nodes)!=0) {
  adj_matrix <- adj_matrix[-disconnected.nodes,-disconnected.nodes]
  groups <- groups[-disconnected.nodes]
  shapes <- shapes[-disconnected.nodes]
  colours <- colours[-disconnected.nodes]
}

n
length(which(colours=='#CC79A7'))
length(which(colours=='#0072B2'))
length(which(colours=='#E69F00'))
length(which(colours=='#CC79A7'))
q

# Create a qgraph with layout options
qgraph_adj_mat <- qgraph(input=adj_matrix,
                         labels=colnames(adj_matrix),
                         groups=groups,
                         DoNotPlot=TRUE,
                         borders=FALSE,
                         palette="colorblind",
                         label.font='sans',
                         posCol="#009E73",# colour of positive edges
                         negCol="#D55E00",      # colour of negative edges
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
                        title=paste0("Full network (delta = ", 
                                     delta,'; limit = ', limit, ')'),
                        collection="ALL+INFECT")

# Extract the first order neighborhoods:
# nei_node_graphs_adj <- make_ego_graph(graph = igraph_adj, order = 1)
# # How many clinical variables are left
# all <- ncol(adj_matrix)
# Plot the first order neighbourhood for each discrete variable
# for (i in 1:all) {
#   total <- names(nei_node_graphs_adj[[i]][1])
#   circle <- layout_in_circle(nei_node_graphs_adj[[i]],
#                              # Remove the first node from the circle
#                              order=total[-which(total==colnames(adj_matrix)[i])])
#   # Adjust the networks when there is a small number of nodes
#   if (length(total)<4) {
#     # Make sure that the vertices are on one line
#     for (coord in 1:nrow(circle)) {
#       if (circle[coord,1]==-1) {
#         circle[coord,2] <- 0
#       }
#     }
#     # Reorder so that a vertical instead of a horizontal network is created
#     circle <- cbind(circle[,2], circle[,1])
#   }
#   # Plot the network
#   plot(nei_node_graphs_adj[[i]], main=colnames(adj_matrix)[i], layout=circle)
# }

# ################################################################################
# 
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
                        title=paste0('First-order neighbourhoods of clinical parameters (delta = ',
                                     delta,'; limit = ', limit, ')'),
                        collection="ALL+INFECT")
# 
# # Extract the first order neighborhoods:
# nei_node_graphs_adj <- make_ego_graph(graph = igraph_adj, order = 1)
# # How many clinical variables are left
# all <- ncol(adj_matrix)
# # Plot the first order neighbourhood for each discrete variable
# for (i in 1:all) {
#   total <- names(nei_node_graphs_adj[[i]][1])
#   circle <- layout_in_circle(nei_node_graphs_adj[[i]],
#                              # Remove the first node from the circle
#                              order=total[-which(total==colnames(adj_matrix)[i])])
#   # Adjust the networks when there is a small number of nodes
#   if (length(total)<4) {
#     # Make sure that the vertices are on one line
#     for (coord in 1:nrow(circle)) {
#       if (circle[coord,1]==-1) {
#         circle[coord,2] <- 0
#       }
#     }
#     # Reorder so that a vertical instead of a horizontal network is created
#     circle <- cbind(circle[,2], circle[,1])
#   }
#   # Plot the network
#   plot(nei_node_graphs_adj[[i]], main=colnames(adj_matrix)[i], layout=circle)
# }

