\name{mni.rows.per.node}
\alias{mni.rows.per.node}
\title{Configure system for parallel processing of dataset}
\description{
  Computes the number of rows that each individual process in a
  paralleized run of mni.vertex.statistics will have to work with.
}
\usage{
mni.rows.per.node(n.rows, n.nodes)
}
\arguments{
  \item{n.rows}{Number of rows in the dataset.}
  \item{n.nodes}{Number of nodes to parallelize across.}
}
\value{
  Returns a matrix with the same number of columns as number of nodes
  given in the n.nodes argument and two columns, the first column
  corresponding to the lower index and the second to the higher index
  into the datatable that is to be used by each process.
}
