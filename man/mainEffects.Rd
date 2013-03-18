\name{mainEffects}
\alias{mainEffects}
\docType{data}
\title{Estimated main effects}
\description{Estimated main effects (single knock-down effects) for target and query genes. Additional overall effects for each phenotype are contained. The dataset contains an annotation of target genes, query genes, and phenotypes.}
\usage{data(mainEffects)}
\format{
  The format is:
List of 4
 $ target: num [1:289, 1:2, 1:2, 1:11, 1:2] -0.31065 -0.32253 -0.08466 -0.00367 -0.60867 ...
  ..- attr(*, "dimnames")=List of 5
  .. ..$ targetGene  : chr [1:289] "TDRD6" "PRDM11" "KDM1B" "INTS12" ...
  .. ..$ targetDesign: chr [1:2] "1" "2"
  .. ..$ queryDesign   : chr [1:2] "1" "2"
  .. ..$ features      : chr [1:11] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
  .. ..$ replicate     : chr [1:2] "1" "2"
 $ query   : num [1:2, 1:20, 1:2, 1:11, 1:2] 0.277 0.265 0.235 0.226 1.165 ...
  ..- attr(*, "dimnames")=List of 5
  .. ..$ targetDesign: chr [1:2] "1" "2"
  .. ..$ queryGene     : chr [1:20] "DPF2" "SMARCA1" "SMARCC1" "SMARCD2" ...
  .. ..$ queryDesign   : chr [1:2] "1" "2"
  .. ..$ features      : chr [1:11] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
  .. ..$ replicate     : chr [1:2] "1" "2"
 $ overall : num [1:2, 1:2, 1:11, 1:2] 0.3685 0.4638 0.3331 0.487 -0.0985 ...
  ..- attr(*, "dimnames")=List of 4
  .. ..$ targetDesign: chr [1:2] "1" "2"
  .. ..$ queryDesign   : chr [1:2] "1" "2"
  .. ..$ features      : chr [1:11] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
  .. ..$ replicate     : chr [1:2] "1" "2"
 $ Anno    :List of 3
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
data(mainEffects, package="HD2013SGI")
print(dim(mainEffects$target))
plot(mainEffects$target[,1,1,1,1])
}
\keyword{datasets}
