\name{mni.read.vertstats.column}
\alias{mni.read.vertstats.column}
\title{Read a column from a vertstats file}
\description{
  Uses a shell command to read a column from a vertstats file.
}
\usage{
column <- mni.read.vertstat.column(filename, column.name=NULL)
}
\arguments{
  \item{filename}{The filename}
  \item{column.name}{The name of the column to read from the vertstats file}
}
\value{
  A numeric array containing the column.
}
\details{
  Uses the shell command "vertstats_extract" to return the specified
  column from the vertstats file. That command thus needs to be on the
  path of the shell from which R was started.
}
\examples{
column <- mni.read.vertstats.column("somefile.vertstats", "tstatistics")
}
\author{Jason Lerch}
\keyword{file}