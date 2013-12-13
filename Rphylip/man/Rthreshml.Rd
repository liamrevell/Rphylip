\name{Rthreshml}
\alias{Rthreshml}
\title{R interface for threshml}
\usage{
Rthreshml(tree, X, types=NULL, path=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{X}{a \emph{data.frame} of continuous valued or discrete character traits with rownames containing species names. Discrete & continuous characters can be supplied in any order. All discrete character traits must be two-state, but can be coded using any convention (i.e., \code{0, 1}, \code{"A","B"}, etc.).}
	\item{types}{character vector containing the types (e.g., \code{"discrete"}, \code{"continuous"}). If \code{types} are not supplied, \code{Rthreshml} will try to figure out which columns via \code{is.numeric}.}
	\item{path}{path to the executable containing threshml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to threshml. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{burnin} burnin generations for the MCMC; \code{nchain} number of chains of the MCMC; \code{ngen} number of generations in each chain; \code{proposal} variance on the proposal distribution for the MCMC; \code{lrtest} logical value indicating whether to conduct a likelihood-ratio test of the hypothesis that some correlations are zero (does not appear to work in the present version); and \code{cleanup} remove PHYLIP input/output files after the analysis is completed (defaults to \code{cleanup = TRUE}).
}
\description{
	This function is an R interface for threshml in the PHYLIP package (Felsenstein 1989, 2013). threshml fits the threshold model of Felsenstein (2005; 2012). Note that threshml is new & not in the currently released version of PHYLIP (as of December 2013). It can be downloaded from its webpage here: \url{http://evolution.gs.washington.edu/phylip/download/threshml/}. If not specifying \code{path}, the executable file for threshml (e.g., threshml.exe in Windows) should be placed in the folder containing all other executable files for PHYLIP (e.g., C:/Program Files/phylip=3.695/exe in Windows).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\value{
	This function returns a list containing the results from threshml.
}
\references{
	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2005) Using the quantitative genetic threshold model for inferences between and within species. \emph{Philosophical Transactions of the Royal Society London B}, 360, 1427–1434.

	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, 179, 145-156.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rcontrast}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
