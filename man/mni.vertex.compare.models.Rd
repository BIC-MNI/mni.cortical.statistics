\name{mni.vertex.compare.models}
\alias{mni.vertex.compare.models}
\title{An ANOVA comparing two models at each vertex}
\description{
  Tests for differences when fitting two models at each vertex.
}
\usage{
mni.vertex.compare.models(glim.matrix, model.one, model.two, vertex.table)
}
\arguments{
  \item{glim.matrix}{The matrix returned by \link{mni.read.glim.file}.}
  \item{model.one}{The first statistical model.}
  \item{model.two}{The second statistical model.}
  \item{vertex.table}{The vertex table returned by \link{mni.build.data.table}.}
}
\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
