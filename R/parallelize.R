mni.rows.per.node <- function(n.rows, n.nodes) {
  # number of rows for each node
  n.per.node <- floor(n.rows / n.nodes)

  # setup the row indices each node will compute over
  node.setup <- matrix(nrow = n.nodes, ncol = 2)
  for (i in 1:n.nodes) {
    node.setup[i,1] <- ((i*n.per.node) - n.per.node) + 1
    node.setup[i,2] <- (i*n.per.node)
  }
  # put any remainders onto the last node
  node.setup[n.nodes,2] = node.setup[n.nodes,2]+n.rows %% n.nodes
  # return the matrix
  return(node.setup)
}

mni.split.data.for.nodes <- function(data.table, node.setup, directory) {
  # setup array of filenames
  filenames <- paste(directory, "/", "node", 1:nrow(node.setup), ".Rdata", sep="")
  # save each segment of the data table in a separate file
  for (i in 1:nrow(node.setup)) {
    node.data <- data.table[node.setup[i,1]:node.setup[i,2],]
    save(node.data, file=filenames[i])
  }
  return(filenames)
}
