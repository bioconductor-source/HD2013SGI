\name{QueryAnnotation}
\alias{QueryAnnotation}
\docType{data}
\title{Annotation of all query genes in the screen}
\description{Annotation of all query genes in the screen.}
\usage{data(QueryAnnotation)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Symbol}}{a character vector}
  }
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(QueryAnnotation, package="HD2013SGI")
print(QueryAnnotation$Symbol)
}
\keyword{datasets}
