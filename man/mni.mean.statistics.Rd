\name{mni.mean.statistics}
\alias{mni.mean.statistics}
\title{Run a linear model on the average of a series of cortical stats files}
\description{
  Uses a GLIM matrix and regresses the average of each file listed in
  the first column against the model.Any formula that is valid for
  \link{lm} is also valid here.
}
\usage{
mni.mean.statistics(glim.matrix, statistics.model, vertex.table=FALSE)
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

\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
