\name{mni.vertex.mixed.model}
\alias{mni.vertex.mixed.model}
\title{Per vertex mixed model regression}
\description{
  Runs a mixed model regression at every vertex.
}
\usage{
vs <- mni.vertex.mixed.model(glim.matrix, fixed.effect, random.effect, vertex.table=FALSE)
}
\arguments{
   \item{glim.matrix}{The design matrix as specified by
     \link{mni.read.glim.file}}
   \item{fixed.effect}{The fixed effect part of the formula}
   \item{random.effect}{The random effects part of the formula}
   \item{vertex.table}{The data table as created by \link{mni.build.data.table}}
}
\value{
  A list with the same number of rows as vertices in the data table
  containing the value, standard error, and t-statistic.
}
\details{
  Runs a linear mixed effect model at each vertex. The two formula
  specifications are the same as in \link{lme}, so that is a good place
  to get some more in depth guidance.
}
\seealso{
  \link{mni.vertex.statistics}
  \link{lme}
}
\examples{
}
\author{Jason Lerch}
\keyword{regression}