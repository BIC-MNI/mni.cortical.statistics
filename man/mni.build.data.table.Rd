\name{mni.build.data.table}
\alias{mni.build.data.table}
\title{Load all vertices from a glim matrix}
\description{
  Loads all the vertex information from disk. This function is called
  automatically in \link{mni.mean.statistics} and
  \link{mni.vertex.statistics}, but can be used separately to speed up
  computation for multiple statistics runs or for the plotting of
  individual vertices.
}
\usage{
mni.build.data.table(glim.matrix)
}
\arguments{
  \item{glim.matrix}{The matrix returned by \link{mni.read.glim.file}.}
}
\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
