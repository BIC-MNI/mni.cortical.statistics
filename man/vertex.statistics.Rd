\name{vertex.statistics}
\alias{vertex.statistics}
\title{Run a linear model at every vertex of a series of cortical stats files}
\description{
  Uses a GLIM matrix and regresses the values at every vertex against
  the linear model.
}
\usage{
vertex.stats <- vertex.statistics(glim.matrix, statistics.model)
}
\arguments{
  \item{glim.matrix}{The matrix to test against. Usually read using
read.glim.file.}
  \item{statistics.model}{The model to test, where y takes the place of
    the independent variable. For example, 'y ~ Intercept + Group'}
}
\seealso{
  \code{\link{read.glim.file}}.
  \code{\link{mean.statistics}}.
  \code{\link{write.vertex.stats}}.
}
