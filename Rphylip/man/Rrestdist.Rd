\name{Rrestdist}
\alias{Rrestdist}
\title{R interface for restdist}
\usage{
Rrestdist(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"rest.data"} containing data for restriction sites.}
	\item{path}{path to the executable containing restdist. If \code{path = NULL}, the R will search several paths for restdist. See details for more information.}
	\item{...}{optional arguments to be passed to restdist. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{method}, the method for calculating the genetic distance - options are \code{"Modified"} (the default), and \code{"Nei/Li"} (Nei & Li 1979); \code{data}, the type of data - either \code{"sites"} or \code{"fragments"}; \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity); \code{kappa} transition:transversion ratio (defaults to \code{kappa = 2.0}); \code{site.length} the restriction site length (defaults to \code{site.length=6}); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the restdist program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/restdist.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for restdist in the PHYLIP package (Felsenstein 2013). restdist can be used to estimate the evolutionary distances between populations based on restriction site or fragment data.
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{
	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Nei, M., Li, W.-H. (1979) Mathematical model for studying genetic variation in terms of restriction endonucleases. \emph{Proceedings of the National Academy of Sciences}, 76, 5269-5273. 
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\seealso{
	\code{\link{Rdnadist}}, \code{\link{Rgendist}}, \code{\link{Rprotdist}}
}
\examples{
\dontrun{
data(restriction.data)
D<-Rrestdist(restriction.data)
}
}
\keyword{phylogenetics}
\keyword{distance method}
