\name{HD2013SGIselectByStability}
\alias{HD2013SGIselectByStability}
\title{Feature selection method}
\description{A function to select features that are most stable across replicated experiments}
\usage{
HD2013SGIselectByStability(subsample, preselect = "count", Rdim = 40, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subsample}{This is the input data. Usually a subsample of the complete screen is enough to select the non-redundant features. \code{subsample} is a list with three elements: \code{D} is a three-dimensional array with dimensions samples x features x replicates. As samples usually 1000 to 5000 experiments are randomly selected. The function needs two replicates.
}
  \item{preselect}{Names of the features that should be preselected, e.g. count is preselected in this screen, because of its biological interpretability and comparability to other viability-based genetic interaction screens.}
  \item{Rdim}{The maximum number of features to be selected.}
  \item{verbose}{If TRUE information about the progress and the quality of the selected features is printed.}
}
\value{
(selected = selected, correlation = correlation, ratioPositive = ratioPositive, correlationAll = correlationAll)
  \item{selected }{The names of the selected features in the order as selected.}
  \item{correlation }{The correlation of the residual features after fitting a linear function on the previously selected features. Correlations are in same order as \code{selected}.}
  \item{ratioPositive }{The fraction of positively correlated features among all candidate features in each step. In same order as \code{selected}.}
  \item{correlationAll }{The correlation of the residual features of all candidate features in each step of the selection process.}
}
\author{Bernd Fischer}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
