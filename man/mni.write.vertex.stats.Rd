\name{mni.write.vertex.stats}
\alias{mni.write.vertex.stats}
\title{Write the output of a vertex statistics run to file}
\description{
  Takes the output from mni.vertex.statistics and writes it to a
  vertstats file using either the new (with headers) or old version of
  the file format.
}
\usage{
mni.write.vertex.stats(vertex.stats, "filename", headers=TRUE,
mean.stats=ms, glim.matrix=gm)
}
\arguments{
  \item{vertex.stats}{The results from mni.vertex.statistics}
  \item{filename}{The filename to write to.}
  \item{headers}{Whether to write out header information. If TRUE, then
    a new style file is created with a hierarchical header; this is the
    default. If false, then a simple file without headers is created.}
  \item{mean.stats}{The results from mni.mean.statistics run on the same
    dataset. This information is written to the header if
    headers=TRUE. mean.stats is NULL by default.}
  \item{glim.matrix}{The glim matrix used with mni.read.glim.file. If
    passed to this argument then the matrix is added to the
    header. Default is NULL}
}
\examples{
# A complete statistical run:
# read the design matrix file
\dontrun{gf <- mni.read.glim.file("filename.glim", header=TRUE)}
# see if the mean is different
\dontrun{ms <- mni.mean.statistics(gf, 'y ~ 1 + V2')}
# look at the results of the mean stats run
\dontrun{summary(ms)}
# run a test at each vertex
\dontrun{vs <- mni.vertex.statistics(gf, 'y ~ 1 + V2')}
# write the results out to file.
\dontrun{mni.write.vertex.stats(vs, "filename.vertstats", headers=TRUE,
mean.stats=ms, glim.matrix=gf)}
}
\keyword{file}
\keyword{regression}
\seealso{
  \code{\link{mni.mean.statistics}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
