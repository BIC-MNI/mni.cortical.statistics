\name{mni.read.glim.file}
\alias{mni.read.glim.file}
\title{Read Linear Model Matrix Text File}
\description{
  Reads a matrix to be used in testing a linear model, where the matrix
  is a simple space separated text file.
}
\usage{
glim.matrix <- mni.read.glim.file(filename, header=FALSE, fill=FALSE)
}
\arguments{
  \item{filename}{The filename of the file containing the matrix.}
  \item{header}{Whether the first row of the matrix contains header
    information.}
  \item{fill}{Whether to allow for missing values or not.}
}
\keyword{file}
\keyword{regression}
\seealso{
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
  \code{\link{mni.build.data.table}}
}
