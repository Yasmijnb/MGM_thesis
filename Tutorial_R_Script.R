myPath <- "../Provided Data/"
Single_cell_Tirosh <- read.csv(paste0(myPath,"GSE72056_melanoma_single_cell_revised_v2.txt"), sep = "\t")
# Remove NAs from the data
#data <- na.omit(data)
# Use a small subset of the data
Single_cell_Tirosh <- Single_cell_Tirosh[1:10000,]

#we restrict our MGM estimation to the top 1000 most abundant genes:
p <- 1000 #number of continuous variables
Annotation <- Single_cell_Tirosh[1:3,]              #encodes the discrete variable "cell type"
                                                    #with 7 different levels
Single_cell_Tirosh <- Single_cell_Tirosh[-c(1:3),]  #original continuous data matrix:
                                                    #rows = genes, columns = samples
gene_names <- make.names(as.character(Single_cell_Tirosh[,1]), unique = TRUE)
Single_cell_Tirosh <- as.matrix(Single_cell_Tirosh[,-1])
rownames(Single_cell_Tirosh) <- gene_names
means <- rowMeans(Single_cell_Tirosh)

#order genes according to RNA abundance:
Single_cell_Tirosh_ordered <- Single_cell_Tirosh[order(means, decreasing = TRUE),]
#take top 1000 genes:
Single_cell_Tirosh_ordered <- Single_cell_Tirosh_ordered[1:p,]
#scale continuous variables:
Single_cell_Tirosh_scaled <- t(scale(t(Single_cell_Tirosh_ordered)))

rownames(Annotation) <- Annotation[,1]
Annotation <- Annotation[,-1]

#discrete variable matrix has to be parameterized as follows:
D <- NULL
for(i in c(0:6)){
  Dtemp <- rep(0,ncol(Annotation))
  Dtemp[as.integer(Annotation[3,])==i] <- 1
  D <- rbind(D,Dtemp)
}

#check if each sample has only one unique level:
all(apply(D,2,sum)==1)

rownames(D) <- c("tumor/unclassif","T cell","B cell","Macro","Endo","CAF","NK")

#each level of a discrete variable is now coded in one individual row, "0" = "not this level", "1" = "this level", the first level will be treated as the baseline:
D[1:7,1:3]

#save matrices as npy-files in Python-script folder (here folder "MGM_algorithm") for uploading in Python:
library(RcppCNPy)
myPath <- "../Try 1/"
setwd(paste0(myPath,"MGM_algorithm/"))
#continuous variable matrix with rows = samples, columns = genes:
npySave("Xsc.npy",t(Single_cell_Tirosh_scaled))
#discrete variable matrix with rows = samples, columns = cell type levels:
npySave("Dsc.npy",t(D))
#vector with numbers of levels for each discrete variable in same order as D:
npySave("levelssc.npy", as.integer(7))
#number of samples
n <- ncol(Single_cell_Tirosh_scaled) 
#number of discrete variables
q <- 1 
#define lambda sequence for penalization according to (Lee and Hastie, 2015):
la_seq <- 5*sqrt(log(p+q)/n)
npySave("lam_seqsc.npy", la_seq) #save lambda sequence

system(paste0("python apply_functions_command_line.py",
              " --iterations=100", #set number of iterations to 1000
              " --oTol=1e-4", #set precision to 1e-4
              " --X_file=Xsc.npy", #continuous variable matrix
              " --D_file=Dsc.npy", #discrete variable matrix
              " --lamseq_file=lam_seqsc.npy", #lambda sequence
              " --levels_file=levelssc.npy", #levels of discrete variables
              " --results_folder=temp"), wait=TRUE) #output folder

setwd(paste0(myPath,"MGM_algorithm/temp/"))

#matrix of edges between continuous - continuous variables:
B0 = npyLoad("B_0.npy")
#matrix of edges between continuous - discrete variables:
Rho0 = npyLoad("Rho_0.npy")
#matrix of edges between discrete - discrete variables:
Phi0 = npyLoad("Phi_0.npy")

#combine single result matrices into one adjacency matrix:
adj_matrix <- rbind(cbind(B0,t(Rho0)), cbind(Rho0,Phi0))
colnames(adj_matrix) <- rownames(adj_matrix) <- c(rownames(Single_cell_Tirosh_scaled),rownames(D))
#color genes and cell type differently
groups <- as.factor(c(rep("gene",p), rep("cell type",7))) 
#shape gene and cell type nodes differently
shapes <- as.factor(c(rep("circle",p), rep("square",7))) 

#remove all nodes which are not connected to any other node in the estimated MGM:
diconnected.nodes <- which(apply(adj_matrix, 1, function(x){all(x==0)}))
adj_matrix <- adj_matrix[-diconnected.nodes,-diconnected.nodes]
groups <- groups[-diconnected.nodes]
shapes <- shapes[-diconnected.nodes]
library(qgraph)
qgraph_adj_mat <- qgraph(adj_matrix,labels = colnames(adj_matrix),
                         groups = groups,
                         DoNotPlot = TRUE,
                         borders=FALSE,
                         palette="colorblind",
                         label.font='sans',
                         posCol=c("royalblue1"),
                         negCol=c("red"),
                         color=c('royalblue1','red1'),
                         shape =as.vector(shapes),
                         fade=TRUE,
                         esize=2)


library(igraph)
#convert to igraph object
igraph_adj <- as.igraph(qgraph_adj_mat, attributes = TRUE) 
#extract the first order neighborhoods:
nei_node_graphs_adj <- make_ego_graph(graph = igraph_adj, order = 1)

#plot the first order neighborhood of cell type "B cell":
plot(nei_node_graphs_adj[[1000]], layout=layout_in_circle(nei_node_graphs_adj[[1000]]))

#plot the first order neighborhood of cell type "T cell":
plot(nei_node_graphs_adj[[999]], layout=layout_in_circle(nei_node_graphs_adj[[999]]))
