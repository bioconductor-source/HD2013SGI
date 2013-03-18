\name{datamatrix}
\alias{datamatrix}
\docType{data}
\title{Phenotypic data after quality control and feature selection}
\description{Phenotypic features of pairwise genetic perturbation experiments \bold{after} selection of non-redundant features and quality control. \code{D} is the 6-dimensional array of experimental measurements. Its dimensions are target genes x target siRNA designs x query genes x query siRNA designs x features x replicates. The array has a dimnames attribute, but there exists a more comprehensive annotation of target genes, query genes, and phenotypes in \code{Anno}.}
\usage{data(datamatrix)}
\format{
  The format is:
List of 2
 $ D   : num [1:289, 1:2, 1:20, 1:2, 1:11, 1:2] -0.05334 -0.20294 -0.10123 0.33203 0.00638 ...
  ..- attr(*, "dimnames")=List of 6
  .. ..$ targetGene  : chr [1:289] "TDRD6" "PRDM11" "KDM1B" "INTS12" ...
  .. ..$ targetDesign: chr [1:2] "1" "2"
  .. ..$ queryGene     : chr [1:20] "DPF2" "SMARCA1" "SMARCC1" "SMARCD2" ...
  .. ..$ queryDesign   : chr [1:2] "1" "2"
  .. ..$ features      : chr [1:11] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
  .. ..$ replicate     : chr [1:2] "1" "2"
 $ Anno:List of 3
  ..$ target :'data.frame':	289 obs. of  4 variables:
  .. ..$ ID    : chr [1:289] "B1" "B2" "B3" "B4" ...
  .. ..$ Symbol: chr [1:289] "TDRD6" "PRDM11" "KDM1B" "INTS12" ...
  .. ..$ Well  : chr [1:289] "B1" "B2" "B3" "B4" ...
  .. ..$ group : chr [1:289] "sample" "sample" "sample" "sample" ...
  ..$ query    :'data.frame':	20 obs. of  2 variables:
  .. ..$ ID    : chr [1:20] "01" "02" "03" "04" ...
  .. ..$ Symbol: chr [1:20] "DPF2" "SMARCA1" "SMARCC1" "SMARCD2" ...
  ..$ phenotype: chr [1:11] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(datamatrix, package="HD2013SGI")
plot(datamatrix$D[,1,1,1,1,1])
}
\keyword{datasets}
