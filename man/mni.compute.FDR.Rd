\name{mni.compute.FDR}
\alias{mni.compute.FDR}
\title{Compute the q value at every vertex}
\description{
  Uses the False Discovery Rate as described by Genovese et. al. to
  compute both the q-values and a t-threshold for a corresponding set of
  t or p values.
}
\usage{
mni.compute.FDR(t.stats=NULL, p.values=NULL, df=Inf, fdr=0.05,
                plot.fdr=FALSE) 
}
\arguments{
  \item{t.stats}{An array of t-statistics.}
  \item{p.values}{An array of p-values.}
  \item{df}{The degrees of freedom.}
  \item{fdr}{The threshold to test for.}
  \item{plot.fdr}{Whether to produce a plot of the FDR procedure.}
}
\value{
  Returns both the threshold as well as a list of q-values corresponding
  to the t or p values.
}
\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
