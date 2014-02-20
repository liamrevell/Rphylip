\name{Rdollop}
\alias{Rdollop}
\title{R interface for dollop}
\usage{
Rdollop(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"phylip.data"} with a binary character; or a matrix. Data should be \code{"0"} and \code{"1"}. The states \code{"P"} (or \code{"B"}), for polymorphism, and \code{"?"} are also permitted.}
	\item{path}{path to the executable containing dollop. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dollop. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the Dollo or polymorphism parsimony score will be computed on a fixed input topology; \code{method} indicating whether Dollo (\code{method="dollo"}, the default) or polymorphism (\code{method="polymorphism"}) should be used; \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{threshold} threshold value for threshold parsimony (defaults to ordinary parsimony); \code{ancestors} vector of ancestral states, which allows for the unordered Dollo parsimony method of Felsenstein (1984); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the dollop program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dollop.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dollop in the PHYLIP package (Felsenstein 2013). dollop can be used for Dollo (Le Quesne 1974; Farris 1977) and polymorphism (Inger 1967; Farris 1978; Felsenstein 1979) parsimony tree inference.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Farris, J.S. (1977) Phylogenetic analysis under Dollo's Law. \emph{Systematic Zoology}, 26, 77-88.

	Farris, J.S. (1978) Inferring phylogenetic trees from chromosome inversion data. \emph{Systematic Zoology}, 27, 275-284.

	Felsenstein, J. (1979) Alternative methods of phylogenetic inference and their interrelationship. \emph{Systematic Zoology}, 28, 49-62.

	Felsenstein, J. (1984) The statistical approach to inferring evolutionary trees and what it tells us about parsimony and compatibility. pp. 169-191 in: \emph{Cladistics: Perspectives in the Reconstruction of Evolutionary History}, edited by T. Duncan and T. F. Stuessy. Columbia University Press, New York. 

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Inger, R.F. (1967) The development of a phylogeny of frogs. \emph{Evolution}, 21, 369-384. 

	Le Quesne, W.J. (1974) The uniquely evolved character concept and its cladistic application. \emph{Systematic Zoology}, 23, 513-517. 
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rpars}}, \code{\link{Rdnapars}}
}
\examples{
\dontrun{
data(primates.bin)
tree<-Rdollop(primates.bin)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
