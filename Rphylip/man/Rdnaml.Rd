\name{Rdnaml}
\alias{Rdnaml}
\title{R interface for dnaml}
\usage{
Rdnaml(X, path="." , ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnaml.}
	\item{...}{optional arguments to be passed to dnaml. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the model will be optimized on a fixed input topology; \code{kappa} transition:transversion ratio (defaults to \code{kappa = 2.0}); \code{bf} vector of base frequencies in alphabetical order (i.e., A, C, G, & T) - if not provided, then defaults to empirical frequencies; \code{rates} vector of rates (defaults to single rate); \code{rate.categories} vector of rate categories corresponding to the order of \code{rates}; \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity); \code{ncat} number of rate categories for the gamma model; \code{inv} proportion of invariant sites for the invariant sites model (defaults to \code{inv = 0}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{speedier} speedier but rougher analysis (defaults to \code{speedier = FALSE}); \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).
}
\description{
	This function is an R interface for dnaml in the PHYLIP package (Felsenstein 2013). dnaml can be used for ML phylogeny estimation from DNA sequences.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{opt.Rdnaml}}, \code{\link{Rcontml}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
