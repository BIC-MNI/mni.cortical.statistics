\name{mni.compute.FDR}
\alias{mni.compute.FDR}
\title{Compute the q value at every vertex}
\description{
  Uses the False Discovery Rate as described by Genovese et. al. to
  compute both the q-values and a t-threshold for a corresponding set of
  t or p values.
}
\usage{
mni.compute.FDR(t.stats=NULL, p.values=NULL, filename=NULL,
                column.name=NULL, df=Inf, fdr=0.05, plot.fdr=FALSE) 
}
\arguments{
  \item{t.stats}{An array of t-statistics.}
  \item{p.values}{An array of p-values.}
  \item{filename}{Quoted filename of previous results}
  \item{column.name}{The column of the file specified by filename to retrieve}
  \item{df}{The degrees of freedom.}
  \item{fdr}{The threshold to test for.}
  \item{plot.fdr}{Whether to produce a plot of the FDR procedure.}
}
\value{
  Returns both the threshold as well as a list of q-values corresponding
  to the t or p values.
}
\details{This function computes the False Discovery Rate as described by
  the NeuroImage paper by Genovese et. al. It can handle three possible
  inputs: an array of p values, and array of t statistics, or a filename
  of vertstats results (as generated by
  \link{mni.write.vertex.stats}). In the latter case, the filename
  should be quoted and accompanied by a column name (specified by the
  column.name argument, also quoted). If you are unsure about the
  possible column names in that file, leave that argument out and
  mni.compute.FDR will report an error but also give a list of the
  possible column names in that file.
}
\keyword{regression}
\seealso{
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
