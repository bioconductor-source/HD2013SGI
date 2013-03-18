\name{HD2013SGI-package}
\alias{HD2013SGI-package}
\alias{HD2013SGI}
\docType{package}
\title{Experimental Data and Analysis of the HCT116 Genetic Interaction Matrix}
\description{Experimental Data and Analysis of the HCT116 Genetic Interaction Matrix. This package contains the data and source code for the paper Laufer, Fischer, Billmann, Huber, Boutros, HD2013SGI, 2013.}
\details{
\tabular{ll}{
Package: \tab HD2013SGI\cr
Type: \tab Package\cr
Version: \tab 0.0.3\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab true\cr
Imports: \tab \cr
Depends: \tab rhdf5, RColorBrewer, gplots, geneplotter, MASS, grid, hwriter, splots, igraph, abind, limma, vcd, LSD\cr
SystemRequirements: \tab GNU make\cr
biocViews: \tab Infrastructure\cr
Built: \tab R 2.15.1; ; 2013-02-14 12:08:09 UTC; unix\cr
}

The interaction matrix can be loaded by \code{data(Interactions, package="HD2013SGI")}. Type \code{?Interactions} to see a documentation of the interaction data.

The vignette of the package can be seen by typing
>library("HD2013SGI")
>vignette("HD2013SGI")
It contains the complete documentation and R-code for the analysis of the data published in the original publication.

All intermediate results are precomputed and can be loaded. the following datasets are available:
\tabular{ll}{
\code{\link{featuresPerWell}} \tab The screen data in screen order\cr
\code{\link{datamatrixfull}} \tab The phenotype data of all pairwise genetic perturbation experiments before quality control and feature selection\cr
\code{\link{QueryAnnotation}} \tab Annotation of all the query genes in the screen\cr
\code{\link{TargetAnnotation}} \tab Annotation of all target genes in the screen\cr
\code{\link{stabilitySelection}} \tab Results from the feature selection step\cr
\code{\link{datamatrix}} \tab The phenotype data of all pairwise genetic perturbation experiments after quality control and feature selection\cr
\code{\link{mainEffects}} \tab \cr estimated main effects (single knock down effects)
\code{\link{nrOfInteractionsPerTarget}} \tab number of interactions per target gene\cr
\code{\link{Interactions}} \tab The genetic interaction data (pi-scores, p-values, and annotation)\cr
}

A number of helper functions are defined in the package and used in the vignette.
\tabular{ll}{
\code{\link{HD2013SGIorderDim}} \tab hclust on one out of three dimensions of an interaction matrix\cr
\code{\link{HD2013SGIHeatmapHuman}} \tab plotting a heatmap of a three dimensional array of pi-scores (target genes x query genes x features)\cr
\code{\link{HD2013SGImaineffects}} \tab estimation main effects\cr
\code{\link{HD2013SGIselectByStability}} \tab feature selection to select features most stable between replicated experiments\cr
}

}
\author{
Bernd Fischer

Maintainer: Bernd Fischer <bernd.fischer@embl.de>
}
\references{
Laufer, Fischer, Billmann, Huber, Boutros, HD2013SGI, 2013.
}
\keyword{ package }
\examples{
data(Interactions, package="HD2013SGI")
}
