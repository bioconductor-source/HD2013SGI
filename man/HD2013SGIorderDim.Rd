\name{HD2013SGIorderDim}
\alias{HD2013SGIorderDim}
\title{hclust on one out of three dimensions of an interaction matrix}
\description{hclust on one out of three dimensions of a three-dimensional array of interaction scores (target genes x query genes x features)}
\usage{
HD2013SGIorderDim(x, i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A three dimensional array to be clustered.}
  \item{i}{The dimension of the array along which the data is clustered.}
}
\value{Returns a cluster hierarchy of class \code{hclst}.}
\author{Bernd Fischer}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
