\name{mni.vertex.power.analysis}
\alias{mni.vertex.power.analysis}
\title{Compute the statistical power at every vertex}
\description{
  Given a standard deviation at every vertex this function computes the
  statistical power both in the form of the n required to recapture a
  change of delta, or the delta required to recapture a change of
  n. Basically a wrapper to \link{power.t.test}.
}
\usage{
mni.vertex.power.analysis(std.file, n=25, alpha=0.005, power=0.995, delta=0.5)
}
\arguments{
  \item{std.file}{A text file with one line per vertex containing the
    standard deviations at each vertex.}
  \item{n}{The number of subjects to compute the power for.}
  \item{alpha}{The threshold at which to compute the power.}
  \item{power}{The power to use in the computation.}
  \item{delta}{The change required in the power computation.}
}
\keyword{regression}
\seealso{
  \code{\link{power.t.test}}.
  \code{\link{mni.read.glim.file}}.
  \code{\link{mni.build.data.table}}.
  \code{\link{mni.vertex.statistics}}.
  \code{\link{mni.write.vertex.stats}}.
}
