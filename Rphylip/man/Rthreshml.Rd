\name{Rthreshml}
\alias{Rthreshml}
\title{R interface for threshml}
\usage{
Rthreshml(tree, X, types=NULL, path=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{X}{a \emph{matrix} of continuous valued traits (in columns) with rownames containing species names. For within-species contrasts analysis, the matrix should contain repeating (identical) row names for conspecifics.}
	\item{types}{character vector containing the types (e.g., \code{"discrete"}, \code{"continuous"}). If \code{types} are not supplied, \code{Rthreshml} will try to figure them out.}
	\item{path}{path to the executable containing contrast. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to threshml. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{burnin} burnin generations for the MCMC; \code{nchain} number of chains of the MCMC; \code{ngen} number of generations in each chain; \code{proposal} variance on the proposal distribution for the MCMC; \code{lrtest} logical value indicating whether to conduct a likelihood-ratio test of the hypothesis that some correlations are zero (does not appear to work in the present version); and \code{cleanup} remove PHYLIP input/output files after the analysis is completed (defaults to \code{cleanup = TRUE}).
}
\description{
	This function is an R interface for threshml in the PHYLIP package (Felsenstein 2013). threshml fits the threshold model of Felsenstein (2005; 2012).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\value{
	This function returns a list containing the results from threshml.
}
\references{

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
