\name{Rdnaml}
\alias{Rdnaml}
\alias{Rdnamlk}
\title{R interfaces for dnaml and dnamlk}
\usage{
Rdnaml(X, path=NULL, ...)
Rdnamlk(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnaml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dnaml or dnamlk. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the model will be optimized on a fixed input topology; \code{kappa} transition:transversion ratio (defaults to \code{kappa = 2.0}); \code{bf} vector of base frequencies in alphabetical order (i.e., A, C, G, & T) - if not provided, then defaults to empirical frequencies; \code{rates} vector of rates (defaults to single rate); \code{rate.categories} vector of rate categories corresponding to the order of \code{rates}; \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity); \code{ncat} number of rate categories for the gamma model; \code{inv} proportion of invariant sites for the invariant sites model (defaults to \code{inv = 0}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{speedier} speedier but rougher analysis (defaults to \code{speedier = FALSE}); \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	Finally \code{clock=TRUE} enforces a molecular clock. The argument \code{clock} is only available for \code{Rdnaml}. If \code{clock=TRUE} then dnamlk is used internally. For \code{Rdnamlk} a molecular clock is assumed, thus \code{Rdnaml(...,clock=TRUE)} and \code{Rdnamlk(...)} are equivalent.

	More information about the dnaml and dnamlk programs in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dnaml.html}, and here \url{http://evolution.genetics.washington.edu/phylip/doc/dnamlk.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnaml in the PHYLIP package (Felsenstein 2013). dnaml can be used for ML phylogeny estimation from DNA sequences (Felsenstein 1981; Felsenstein & Churchill 1996).
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (1981) Evolutionary trees from DNA sequences: A Maximum Likelihood approach. \emph{Journal of Molecular Evolution}, 17, 368-376.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Felsenstein, J., Churchill, G. A. (1996) A Hidden Markov Model approach to variation among sites in rate of evolution. \emph{Molecular Biology and Evolution}, 13, 93-104.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{opt.Rdnaml}}, \code{\link{Rcontml}}, \code{\link{Rproml}}
}
\examples{
\dontrun{
data(primates)
tree<-Rdnaml(primates,kappa=10)
clockTree<-Rdnamlk(primates,kappa=10)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
