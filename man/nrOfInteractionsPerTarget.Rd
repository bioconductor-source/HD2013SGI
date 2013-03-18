\name{nrOfInteractionsPerTarget}
\alias{nrOfInteractionsPerTarget}
\docType{data}
\title{Number of interactions per target gene}
\description{Number of genetic interactions per target gene.}
\usage{data(nrOfInteractionsPerTarget)}
\format{
  The format is:
 int [1:282] 1 3 0 2 0 3 1 1 0 1 ...
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(nrOfInteractionsPerTarget, package="HD2013SGI")
plot(nrOfInteractionsPerTarget)
}
\keyword{datasets}
