\name{stabilitySelection}
\alias{stabilitySelection}
\docType{data}
\title{Results from the feature selection method}
\description{Results from the feature selection method.}
\usage{data(stabilitySelection)}
\format{
  The format is:
List of 4
 $ selected      : chr [1:25] "count" "cell.act.m.majoraxis" "nuc.nuc.b.q001" "nuc.0.m.eccentricity" ...
 $ correlation   : num [1:25] 0.917 0.972 0.938 0.928 0.896 ...
 $ ratioPositive : num [1:25] 1 1 0.947 0.942 0.937 ...
 $ correlationAll: List of 25
  ..$ : Named num [1:227] 0.917 0.884 0.93 0.897 0.882 ...
  .. ..- attr(*, "names")= chr [1:227] "count" "nuc.0.m.majoraxis" "nuc.0.m.eccentricity" "nuc.0.s.area" ...
  ..$ : Named num [1:226] 0.884 0.931 0.899 0.887 0.9 ...
  .. ..- attr(*, "names")= chr [1:226] "nuc.0.m.majoraxis" "nuc.0.m.eccentricity" "nuc.0.s.area" "nuc.0.s.perimeter" ...
  ..$ : Named num [1:225] 0.884 0.934 0.884 0.882 0.883 ...
  .. ..- attr(*, "names")= chr [1:225] "nuc.0.m.majoraxis" "nuc.0.m.eccentricity" "nuc.0.s.area" "nuc.0.s.perimeter" ...
  ...
}
\details{\code{selected} is a vector of the selected feature namess. \code{correlation} are the Pearson correlation coefficients of the residual features. \code{ratioPositive} is the fraction of positively correlated features among all candidate features for selection. \code{correlationAll} contains a vector of correlations of the residual features of all candidate features for each step in the selection process.
}
\references{Laufer, Fischer et al., 2013}
\seealso{
\code{\link{HD2013SGI}}
}
\examples{
data(stabilitySelection, package="HD2013SGI")
barplot(stabilitySelection$correlation,
        names.arg=stabilitySelection$selected,las=2)
barplot(stabilitySelection$ratioPositive-0.5, offset=0.5,
        names.arg=stabilitySelection$selected,las=2)
}
\keyword{datasets}
