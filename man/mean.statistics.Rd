\name{mean.statistics}
\alias{mean.statistics}
\title{Run a linear model on the average of a series of cortical stats files}
\description{
  Uses a GLIM matrix and regresses the average of each file listed in
  the first column against the model. 
}
\usage{
mean.stats <- mean.statistics(glim.matrix, statistics.model)
}
\arguments{
  \item{glim.matrix}{The matrix to test against. Usually read using
read.glim.file.}
  \item{statistics.model}{The model to test, where y takes the place of
    the independent variable. For example, 'y ~ Intercept + Group'}
}
\seealso{
  \code{\link{read.glim.file}}.
  \code{\link{vertex.statistics}}.
  \code{\link{write.vertex.stats}}.
}
