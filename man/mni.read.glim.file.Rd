\name{mni.read.glim.file}
\alias{mni.read.glim.file}
\title{Read Linear Model Matrix Text File}
\description{
  Reads a matrix to be used in testing a linear model, where the matrix
  is a simple space separated text file.
}
\usage{
glim.matrix <- mni.read.glim.file(filename, header=FALSE, fill=FALSE,
file.type="space") 
}
\arguments{
  \item{filename}{The filename of the file containing the matrix.}
  \item{header}{Whether the first row of the matrix contains header
    information.}
  \item{fill}{Whether to allow for missing values or not.}
  \item{file.type}{Either "space" or "csv".}
}
\details{
  There are two key ways to control a run of cortical analysis:
  specification of the input data and specification of the formula. The
  input data, read by this function, is prepared as a simple text file
  with one line per subject, and one column per variable. The first
  column is a filename containing the vertstats information, all
  subsequent columns are independent variables. Names of factors are
  allowed, as are missing values if the fill options is set to TRUE. The
  first line can also optionally be a header, set the header option to
  TRUE if that is the case. The different columns are separated by
  either a spaces or commas - space is the default, specify
  file.type="csv" if the file has comma separated values.
}
\value{If there is a header, the returned data
  frame will take the name assigned by the header. If there is no
  header, the column names will be V1, V2, V3 and so forth.}
\keyword{file}
\keyword{regression}
\seealso{
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
  \code{\link{mni.build.data.table}}
}
