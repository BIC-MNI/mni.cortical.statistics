\name{read.glim.file}
\alias{read.glim.file}
\title{Read Linear Model Matrix Text File}
\description{
  Reads a matrix to be used in testing a linear model, where the matrix
  is a simple space separated text file.
}
\usage{
glim.matrix <- read.glim.file(filename, header=FALSE)
}
\arguments{
  \item{filename}{The filename of the file containing the matrix.}
  \item{header}{Whether the first row of the matrix contains header
    information.}
}
\seealso{
  \code{\link{mean.statistics}}.
  \code{\link{vertex.statistics}}.
  \code{\link{write.vertex.stats}}.
}
