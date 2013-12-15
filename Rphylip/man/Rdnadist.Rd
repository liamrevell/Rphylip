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

	More information about the dnadist program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dnadist.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnadist in the PHYLIP package (Felsenstein 2013). dnadist can be used to estimate the evolutionary distances between DNA sequences under various models.
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{

	Barry, D., Hartigan, J.A. (1987) Statistical analysis of hominoid molecular evolution. \emph{Statistical Science}, 2, 191-200.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Felsenstein, J., Churchill, G. A. (1996) A Hidden Markov Model approach to variation among sites in rate of evolution. \emph{Molecular Biology and Evolution}, 13, 93-104.

	Jukes, T.H., Cantor, C.R. (1969) Evolution of protein molecules. pp. 21-132 in \emph{Mammalian Protein Metabolism Vol. III}, ed. M.N. Munro. Academic Press, New York.

	Kimura, M. (1980) A simple model for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. \emph{Journal of Molecular Evolution}, 16, 111-120.

	Kishino, H., Hasegawa, M. (1989) Evaluation of the maximum likelihood estimate of teh evolutionary tree topology from DNA sequence data, and the branching order in Hominoidea. \emph{Journal of Molecular Evolutioon}, 29, 170-179.

	Lake, J.A. (1994) Reconstructing evolutionary trees from DNA and protein sequences: Paralinear distances. \emph{Proceedings of the National Academy of Sciences}, 91, 1455-1459.

	Lockhart, P.J., Steel, M.A., Hendy, M.D., Penny, D. (1994) Recovering evolutionary trees under a more realistic model of sequence evolution. \emph{Molecular Biology and Evolution}, 11, 605-612.

	Steel, M.A. (1994) Recovering a tree from the Markov leaf colourations it generates under a Markov model. \emph{Applied Mathematics Letters}, 7, 19-23.
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
