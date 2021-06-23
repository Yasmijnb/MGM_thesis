# Make an edge table for the cytokines
# Don't remove disconnected nodes before this

# Get the names of all cytokines
cytokines <- colnames(cytokine_data)
cytokines <- cytokines[-(1:13)]

# Make some adjustments
for (odd in c('S100A9', 'ICAM1', 'S100A8')) {
  cytokines[which(cytokines==odd)] <- paste0(odd, '.x')
}

# Make an empty data frame
cyto.edge <- data.frame()

# Fill the data frame as an edge list
for (cyto in 1:length(cytokines)) {
  edges <- NULL
  location <- which(colnames(adj_matrix)==cytokines[cyto])
  for (entry in 1:ncol(adj_matrix)) {
    if (adj_matrix[location,entry] != 0) {
      edges <- c(edges, colnames(adj_matrix)[entry])
      cyto.edge[cyto, length(edges)] <- colnames(adj_matrix)[entry]
    }
  }
}
# Add rownames
rownames(cyto.edge) <- cytokines

View(cyto.edge)

# Save the table as a csv file
# write.csv(cyto.edge, '../Edge list of cytokines (bacterial data).csv')
