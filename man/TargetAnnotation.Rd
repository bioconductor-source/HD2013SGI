\name{TargetAnnotation}
\alias{TargetAnnotation}
\docType{data}
\title{Annotation of all target genes in the screen}
\description{Annotation of the  target genes on one target plate. It includes an ENSEMBL gene identifier, the HUGO name, the position on the plate (well), and the group of the target siRNA (sample or control).}
\usage{data(TargetAnnotation)}
\format{
  A data frame with 345 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Symbol}}{a character vector}
    \item{\code{Well}}{a character vector}
    \item{\code{group}}{a character vector}
  }
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(TargetAnnotation, package="HD2013SGI")
print(TargetAnnotation$Symbol)
}
\keyword{datasets}
