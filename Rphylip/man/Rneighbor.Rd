\name{Rneighbor}
\alias{Rneighbor}
\title{R interface for neighbor}
\usage{
Rneighbor(D, path=NULL , ...)
}
\arguments{
	\item{D}{a distance matrix as an object of class \code{"matrix"} or \code{"dist"}. If a matrix, then \code{D} should be symmetrical with a diagonal of zeros.}
	\item{path}{path to the executable containing neighbor. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to neighbor. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{method} - can be \code{"NJ"} or \code{"nj"} (for neighbor-joining), or \code{"UPGMA"} or \code{"UPGMA"} (for UPGMA); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired (only works for \code{method = "NJ"}, UPGMA trees are already rooted); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for neighbor in the PHYLIP package (Felsenstein 2013). neighbor can be used for neighbor-joining and UPGMA phylogeny inference.
}
\value{
	This function returns an object of class \code{"phylo"} that is the NJ or UPGMA tree.
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnaml}}
}
\examples{
\dontrun{
data(primates)
D<-dist.dna(data(primates),model="JC")
tree<-Rneighbor(D)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{distance method}
