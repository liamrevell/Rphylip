\name{Rdnacomp}
\alias{Rdnacomp}
\title{R interface for dnacomp}
\usage{
Rdnacomp(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnacomp. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dnacomp. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the number of compatible sites will be computed on a fixed input topology; \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the dnacomp program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dnacomp.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnacomp in the PHYLIP package (Felsenstein 2013). dnacomp can be used for phylogeny inference from DNA sequences using the compatibility method (Le Quesne 1969; Fitch 1975).
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Fitch, W.M. (1975) Toward finding the tree of maximum parsimony. pp. 189-230 in \emph{Proceedings of the Eighth International Conference on Numerical Taxonomy}, ed. G. F. Estabrook. W. H. Freeman, San Francisco.

	Le Quesne, W.J. (1969) A method of selection of characters in numerical taxonomy. \emph{Systematic Zoology}, 18, 201-205. 
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnapars}}, \code{\link{Rdnapenny}}
}
\examples{
\dontrun{
data(primates)
tree<-Rdnacomp(primates)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
