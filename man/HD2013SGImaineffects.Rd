\name{HD2013SGImaineffects}
\alias{HD2013SGImaineffects}
\title{Estimation of main effects}
\description{A function to estimate main effects (single knock-down effects) in genetic interaction screens.}
\usage{
HD2013SGImaineffects(x, TP, TargetNeg, QueryNeg, eps = 1e-04, maxiter = 100, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Two dimensional array.}
  \item{TP}{Assignment of target genes to target plates. Used to compute target main effects for each target plate separately.}
  \item{TargetNeg}{Negative controls within the set of target genes.}
  \item{QueryNeg}{Negative controls within the set of query genes.}
  \item{eps}{real number greater than 0. A tolerance for convergence.}
  \item{maxiter}{the maximum number of iterations}
  \item{na.rm}{logical. Should missing values be removed?}
}
\value{
  \item{neg }{Effect of the negative control.}
  \item{targetMainEffect }{target main effects}
  \item{queryMainEffect }{query main effects}
  \item{pi }{Pairwise interaction scores (pi-scores)}
}
\author{changes applied by Bernd Fischer to the implementation of R stats function\code{\link{medpolish}}}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
