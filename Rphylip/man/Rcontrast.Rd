\name{Rcontrast}
\alias{Rcontrast}
\title{R interface for contrast}
\usage{
Rcontrast(tree, X, path=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{X}{a \emph{matrix} of continuous valued traits (in columns) with rownames containing species names. For one trait, \code{X} can be a matrix with one column or a vector with \code{names(X)} containing species names matching \code{tree$tip.label}. For within-species contrasts analysis, the matrix should contain repeating (identical) row names for conspecifics.}
	\item{path}{path to the executable containing contrast. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to contrast. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); and \code{cleanup} remove PHYLIP input/output files after the analysis is completed (defaults to \code{cleanup = TRUE}).
}
\description{
	This function is an R interface for contrast in the PHYLIP package (Felsenstein 2013). contrast can be used to perform the among species phylogenetically independent contrasts method of Felsenstein (1985) and the within & among species method of Felsenstein (2008).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\value{
	If \code{X} contains one observation per species (say, the species mean), then \code{Rcontrast} returns a list with the following components: \code{Contrasts}, a matrix with all phylogenetically independent contrats; \code{Covariance_matrix}, a matrix containing the evolutionary variances (on diagonals) and covariances; \code{Regressions}, a matrix containing the pair-wise bivariate regression coefficients (columns on rows); \code{Correlations}, a correlation matrix of contrasts.

	If \code{X} contains more than one sample per species, then \code{Rcontrast} returns a list with the following elements: \code{VarA}, the estimated among-species variance-covariance matrix; \code{VarE}, the estimated within-species (i.e., 'environmental') variance-covariance matrix; \code{VarA.Regression}, a matrix containing the pair-wise bivariate among-species regression coefficients (columns on rows); \code{VarA.Correlations}, a matrix with the among-species evolutionary correlations; \code{VarE.Regressions}, the pair-wise bivariate within-species regression coefficients; \code{VarE.Correlations}, the within-species correlations; \code{nonVa.VarE}, \code{nonVa.VarE.Regressions}, and \code{nonVa.VarA.Correlations}, estimates obtained when \code{VarA} is not included in the model; \code{logLik} and \code{nonVa.logLik}, log-likelihood when \code{VarA} is included (or not) in the model; \code{k} and \code{nonVa.k} the number of parameters estimated in each model; and \code{P} the p-value of a likelihood-ratio test of \code{VarA}, in which \code{df = k - nonVa.k}. 
}
\references{
	Felsenstein, J. (1985) Phylogenies and the comparative method. American Naturalist, 125, 1-15.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2008) Comparative methods with sampling error and within-species variation: Contrasts revisited and revised. American Naturalist, 171, 713-725.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{pic}}, \code{Rcontml}, \code{Rthreshml}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{continuous characters}
\keyword{maximum likelihood}
