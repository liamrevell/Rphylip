\name{Rclique}
\alias{Rclique}
\title{R interface for clique}
\usage{
Rclique(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"phylip.data"} with a set of binary characters; or a matrix.}
	\item{path}{path to the executable containing clique. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to clique. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{ancestral} vector of ancestral states; \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{minimum.clique}, if specified, will cause \code{Rclique} to print out all cliques (and corresponding trees) greater than \code{minimum.clique}; \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the clique program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/clique.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for clique in the PHYLIP package (Felsenstein 2013). clique can be used for phylogeny inference using the compatibility method (Le Quesne 1969; Estabrook et al. 1976).
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Estabrook, G. F., Johnson,  C. S. Jr., McMorris, F. R. (1976). A mathematical foundation for the analysis of character compatibility. \emph{Mathematical Biosciences}, 23, 181-187. 

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Le Quesne, W.J. (1969) A method of selection of characters in numerical taxonomy. \emph{Systematic Zoology}, 18, 201-205. 
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnacomp}}, \code{\link{Rmix}}
}
\examples{
\dontrun{
data(primates.bin)
tree<-Rclique(primates.bin)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
