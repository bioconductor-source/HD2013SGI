\name{datamatrixfull}
\alias{datamatrixfull}
\docType{data}
\title{Phenotypic data before quality control and feature selection}
\description{Phenotypic features of pairwise genetic perturbation experiments \bold{before} selection of non-redundant features and quality control. \code{D} is the 6-dimensional array of experimental measurements. Its dimensions are target genes x target siRNA designs x query genes x query siRNA designs x features x replicates. The array has a dimnames attribute.}
\usage{data(datamatrixfull)}
\format{
  The format is:
List of 1
 $ D: num [1:345, 1:2, 1:20, 1:2, 1:353, 1:2] 2686 2573 2650 3000 2733 ...
  ..- attr(*, "dimnames")=List of 6
  .. ..$ targetGene  : chr [1:345] "B1" "B2" "B3" "B4" ...
  .. ..$ targetDesign: chr [1:2] "1" "2"
  .. ..$ queryGene     : chr [1:20] "01" "02" "03" "04" ...
  .. ..$ queryDesign   : chr [1:2] "1" "2"
  .. ..$ features      : chr [1:353] "count" "nuc.0.m.cx" "nuc.0.m.cy" "nuc.0.m.majoraxis" ...
  .. ..$ replicate     : chr [1:2] "1" "2"
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(datamatrixfull, package="HD2013SGI")
plot(datamatrixfull$D[,1,1,1,1,1])
}
\keyword{datasets}
