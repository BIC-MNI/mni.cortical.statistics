\name{mni.vertex.statistics}
\alias{mni.vertex.statistics}
\title{Run a linear model at every vertex of a series of cortical stats files}
\description{
  Uses a GLIM matrix and regresses the values at every vertex against
  the linear model. Any formula that is valid for \link{lm} is also
  valid here. 
}
\usage{
mni.vertex.statistics(glim.matrix, statistics.model, vertex.table=FALSE)
}
\arguments{
  \item{glim.matrix}{The matrix to test against. Usually read using
    mni.read.glim.file.}
  \item{statistics.model}{The model to test, where y takes the place of
    the independent variable. For example, 'y ~ Intercept + Group'}
  \item{vertex.table}{A matrix containing the vertex data on which the
    analysis will be run. This matrix is best generated through
    \link{mni.build.data.table}. Alternately, the argument can be left
    as FALSE, in which case the vertex table will be built from the glim
    matrix specified as the first argument.}
}
\value{
  Returns a list, each row corresponding to a vertex. There are columns
  corresponding to all the major statistical terms, such as the adjusted
  R squared, the fstatistics, etc. This list is usually fed into
  \link{mni.write.vertex.stats} for output to file.
}
\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
