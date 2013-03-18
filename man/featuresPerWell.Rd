\name{featuresPerWell}
\alias{featuresPerWell}
\docType{data}
\title{Original phenotypic measurements in screen order}
\description{Original phenotypic features in screen order. \code{Anno} contains the annotation for each experiment including the plate name, row, col, and field. \code{data} is a \code{data.frame} with a column for each phenotypic feature and rows as much as there are experiments in the screen.}
\usage{data(featuresPerWell)}
\format{
  The format is:
List of 2
 $ Anno:'data.frame':	231840 obs. of  4 variables:
  ..$ plate: chr [1:231840] "001CIQ01IRI" "001CIQ01IRI" "001CIQ01IRI" "001CIQ01IRI" ...
  ..$ row  : chr [1:231840] "B" "B" "B" "B" ...
  ..$ col  : chr [1:231840] "1" "1" "1" "1" ...
  ..$ field: chr [1:231840] "1" "2" "3" "4" ...
 $ data: num [1:231840, 1:353] 2780 3120 2242 2603 2170 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:353] "count" "nuc.0.m.cx" "nuc.0.m.cy" "nuc.0.m.majoraxis" ...
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(featuresPerWell, package="HD2013SGI")
plot(log2(featuresPerWell$data[,1]),pch=".")
}
\keyword{datasets}
