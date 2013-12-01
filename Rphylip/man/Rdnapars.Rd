\name{Rdnapars}
\alias{Rdnapars}
\title{R interface for dnapars}
\usage{
Rdnapars(X, path="." , ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnapars.}
	\item{...}{optional arguments to be passed to dnapars. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the parsimony score will be computed on a fixed input topology; \code{thorough} logical value indicating whether to conduct a more thorough search (defaults to \code{thorough=TRUE}); \code{nsave} number of trees to save if multiple equally parsimonious trees are found (defaults to \code{nsave=10000}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{threshold} threshold value for threshold parsimony (defaults to ordinary parsimony); \code{transversion} logical value indicating whether to use transversion parsimony (defaults to \code{transversion=FALSE}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).   
}
\description{
	This function is an R interface for dnapars in the PHYLIP package (Felsenstein 2013). dnapars can be used for ML phylogeny estimation from DNA sequences.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnaml}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
