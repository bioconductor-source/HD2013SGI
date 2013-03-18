\name{HD2013SGIHeatmapHuman}
\alias{HD2013SGIHeatmapHuman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting heatmaps of genetic interaction scores}
\description{Plotting heatmaps of three-dimensional arrays of interaction scores. Two dimensions of the array will be flattened.}
\usage{
HD2013SGIHeatmapHuman(x, cuts, col, colnames = TRUE, rownames = FALSE, mrow = 10, mcol = 10, cexrow = 1, cexcol = 1, border = 0.1, space = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A three dimensional array to be plotted as a heatmap.}
  \item{cuts}{cuts on the values of x for color coding. \code{length(cuts)} has be one larger than \code{length(col)}.}
  \item{col}{Values of \code{x} are mapped on color definitions as defined in \code{col} using the \code{cuts} argument.}
  \item{colnames}{Gene names for columns of the matrix.}
  \item{rownames}{Gene names for rows of the matrix.}
  \item{mrow}{row margin for printing gene names..}
  \item{mcol}{column margin for printing gene names.}
  \item{cexrow}{cex for the row names.}
  \item{cexcol}{cex for the column names.}
  \item{border}{line width of the border.}
  \item{space}{spacing between elements of the third array dimension of \code{x} after flattening.}
}
\value{Nothing is returned.}
\author{Bernd Fischer}
\references{Laufer, Fischer et al., 2013}

\seealso{
\code{\link{HD2013SGI}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
