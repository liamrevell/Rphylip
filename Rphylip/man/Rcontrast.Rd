\name{Rcontrast}
\alias{Rcontrast}
\title{R interface for contrast}
\usage{
Rcontrast(tree, X, path=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{X}{a \emph{matrix} of continuous valued traits (in columns) with rownames containing species names. For within-species contrasts analysis, the matrix should contain repeating (identical) row names for conspecifics.}
	\item{path}{path to the executable containing contrast. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file. (Currently works only for \code{.Platform$OS.type = "windows"}.)}
	\item{...}{optional arguments to be passed to contrast. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); and \code{cleanup} remove PHYLIP input/output files after the analysis is completed (defaults to \code{cleanup = TRUE}).
}
\description{
	This function is an R interface for contrast in the PHYLIP package (Felsenstein 2013). contrast can be used to perform the among species phylogenetically independent contrasts method of Felsenstein (1985) and the within & among species method of Felsenstein (2008).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (1985) Phylogenies and the comparative method. American Naturalist, 125, 1-15.

	Felsenstein, J. (2008) Comparative methods with sampling error and within-species variation: Contrasts revisited and revised. American Naturalist, 171, 713-725.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{pic}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
