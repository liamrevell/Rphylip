\name{Rdnadist}
\alias{Rdnadist}
\title{R interfaces for dnadist}
\usage{
Rdnadist(X, method=c("F84","K80","JC","LogDet"), path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{method}{method for calculating the distances. Can be \code{"F84"} (Kishino & Hasegawa 1989; Felsenstein & Churchill 1996), \code{"K80"} (Kimura 1980), \code{"JC"} (Jukes & Cantor 1969), or \code{"LogDet"} (Barry & Hartigan 1987; Lake 1994; Steel 1994; Lockhart et. al. 1994). Also \code{method="similarity"} computes the sequence similarity among the rows of \code{X}.}
	\item{path}{path to the executable containing dnadist. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dnadist. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity); \code{kappa} transition:transversion ratio (defaults to \code{kappa = 2.0}); \code{rates} vector of rates (defaults to single rate); \code{rate.categories} vector of rate categories corresponding to the order of \code{rates}; \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{bf} vector of base frequencies in alphabetical order (i.e., A, C, G, & T) - if not provided, then defaults to empirical frequencies; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnadist in the PHYLIP package (Felsenstein 2013). dnadist can be used to estimate the evolutionary distances between DNA sequences under various models.
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rneighbor}}
}
\examples{
\dontrun{
data(primates)
D<-Rdnadist(primates,kappa=10)
tree<-Rneighbor(D)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
\keyword{distance method}
