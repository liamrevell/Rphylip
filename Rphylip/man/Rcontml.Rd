\name{Rcontml}
\alias{Rcontml}
\title{R interface for contml}
\usage{
Rcontml(X, path=NULL, ...)
}
\arguments{
	\item{X}{either (a) a \emph{matrix} of continuous valued traits (in columns) with rownames containing species names; or (b) a list of matrices in which each row contains the relative frequency of alleles at a locus for a species. In the latter case the rownames of each matrix in the list should contain the species names.}
	\item{path}{path to the executable containing contml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to contml. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the model will be optimized on a fixed input topology; \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input/output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for contml in the PHYLIP package (Felsenstein 2013). contml can be used for ML phylogeny estimation from gene frequency data or continuous characters. The continuous characters should be rotated so as to be uncorrelated before analysis.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (1981) Maximum likelihood estimation of evolutionary trees from continuous characters. \emph{American Journal of Human Genetics}, 25, 471-492.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.
	
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnaml}}, \code{\link{Rproml}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
\keyword{continuous characters}
