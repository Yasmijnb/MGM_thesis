################################################################################

# Load all packages
library(RcppCNPy)
library(qgraph)
library(igraph)
library(RCy3)
library(stringr)
library(missForest)

################################################################################
## LOAD THE DATA

# Load the data
data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/Merged/Prepped_INFECT.csv",
                 check.names = FALSE)
rownames(data) <- data$Row.names
data <- data[,-c(1:3)]

# Keep only columns which have less than a limit missing data
to_keep <- NULL
limit <- 10 # percentage of missing data that is allowed per variable
for (i in 1:ncol(data)){
  percentage <- 100*sum(is.na(data[,i]))/nrow(data)
  # print(paste(i,':',round(percentage),colnames(data)[i]))
  if (percentage<=limit) {
    to_keep <- c(to_keep, i)
  }
}
if (limit < 100) {
data <- data[,to_keep]
}

################################################################################

# Print some patient characteristics
n.male <- length(which(data$`Sex (0=female, 1=male)` == 1))
p.male <- n.male*100/nrow(data)
n.ss <- length(which(data$`Septic shock Baseline, 1=yes, 0=no` == 1))
p.ss <- n.ss*100/nrow(data)
# I don't know if this is correct
n.amp <- length(which(data$`Amputation of Limb (Time: from diagnosis to ICU day 7)` == 'yes'))
p.amp <- n.amp*100/nrow(data)
n.death <- length(which(data$`Discharged To (1=nursing home, 2=rehab. facil., 3=home, 4=death)(if not in hospital day 90)` == 4))
p.death <- n.death*100/nrow(data)

pat.char <- rbind(c(n.male, p.male),
                  c(n.ss, p.ss),
                  c(n.amp, p.amp),
                  c(n.death, p.death))
colnames(pat.char) <- c('N', '%')
rownames(pat.char) <- c('Sex (male)', 'Septic shock', 'Amputation', 'Death')
pat.char
print(paste(round(mean(data$`Age (years)`),2), '±', round(sd(data$`Age (years)`),2)))

################################################################################
## SPLIT THE DATA INTO DISCRETE AND CONTINUOUS

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
for (column in c('korrektion','HOME','Alcohol Consumption',"foot_amputated",
                 'Mechanical ventilation Baseline',"toe_amputated",
                 "Need For reconstructive Surgery ((if not in hospital day 90)",
                 "finger_amputated","lower_arm_amputated","upper_arm_amputated",
                 "upper_leg_amputated","hand_amputated",
                 "lower_leg_amputated","penis_amputated")) {
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
}
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
# If the variable has only 1 level, remove it from the data
if (is.null(remove_level)==FALSE) {
  Annotation <- Annotation[,-remove_level]
  levels <- levels[-remove_level]
}

# Impute the missing data
data_imp <- missForest(xmis = as.matrix(data), verbose = TRUE, ntree = 100)
data_imp$OOBerror
data <- as.data.frame(data_imp$ximp)

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
  name <- str_remove_all(string = colnames(Annotation)[colnum], pattern = '\\((\\d+\\=[\\w\\s\\-]+[\\,\\s\\)]+)+')
  # If the variable has names as levels
  # if (is.null(levels(Annotation[,colnum]))==FALSE) {
  if (class(Annotation[,colnum])=='character') {
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
    }
  } else
  # If the variable has numbers as levels
    for(i in 0:(levels[colnum]-1)){
      # Start at the minimum level
      mini <- min(Annotation[,colnum])
      i <- i+mini
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
        if ('no' %in% splitted) {
          splitted <- c(splitted, 2, 'unknown')
        }
        level<-splitted[which(splitted==i)+1]
        if (level == 'no' || level == 'yes') {
          colnames_D <- c(colnames_D,name)
        }
        if (level != 'no' && level != 'yes') {
          colnames_D <- c(colnames_D,paste(name,'=',level))
        }
      } 
      else { # if (levels[colnum]!=2)
        colnames_D <- c(colnames_D,paste(colnames(Annotation)[colnum],'=',i))}
  }
}
# Transpose D so it looks the same as the original dataset 
# (samples as rows, variables as columns)
D <- as.data.frame(t(D))
# Add names to the columns
colnames(D) <- colnames_D
colnames(D)[(ncol(D)-10):ncol(D)]

################################################################################
## RUN THE SIMULATION

# Save matrices as npy-files in Python-script folder (here: "MGM_algorithm")
myPath <- "~/School/WUR/SSB-80336 - Thesis/Try 1/MGM_algorithm/"
setwd(myPath)
# Continuous variable matrix with rows = samples, columns = genes:
npySave("Xsc.npy", as.matrix(data))
# Discrete variable matrix with rows = samples, columns = discrete  levels:
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
delta <- 10
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
## MAKE THE NETWORK

# Combine single result matrices into one adjacency matrix:
adj_matrix <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(adj_matrix) <- rownames(adj_matrix) <- c(colnames(data),colnames(D))

# Make different groups
groups <- as.vector(c(rep("continuous",ncol(data)), rep("discrete",ncol(D))))
# Shape discrete variables nodes differently
shapes <- as.vector(c(rep("circle",ncol(data)), rep("circle",ncol(D))))
# Colour discrete variables nodes differently
colours <- as.vector(c(rep("#56B4E9",ncol(data)), rep("#56B4E9",ncol(D))))

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

# Connect to cytoscape (Make sure cytoscape is opened)s
cytoscapePing()
# Create the network
createNetworkFromIgraph(igraph = igraph_adj,
                        title=paste0("Full network (lambda x ", delta,
                                     '; limit = ',limit,')'),
                        collection="INFECT")
