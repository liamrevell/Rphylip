\name{Rdnainvar}
\alias{Rdnainvar}
\title{R interface for dnainvar}
\usage{
Rdnainvar(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}. \code{X} cannot contain more than four species.}
	\item{path}{path to the executable containing dnainvar. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dnainvar. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the dnainvar program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dnainvar.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnainvar in the PHYLIP package (Felsenstein 2013). dnainvar reads DNA sequences and computest Lake's invariants (Cavender & Felsenstein 1987; Lake 1987).
}
\value{
	At present, this function merely prints the results of dnainvar to screen.
}
\references{
	 Cavender, J.A., Felsenstein, J. (1987) Invariants of phylogenies in a simple case with discrete states. \emph{Journal of Classification}, 4, 57-71. 

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	 Lake, J.A. (1987) A rate-independent technique for analysis of nucleic acid sequences: evolutionary parsimony. \emph{Molecular Biology and Evolution}, 4, 167-191. 
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\seealso{
	\code{\link{Rdnapars}}
}
\examples{
\dontrun{
data(primates)
primates<-primates[sample(nrow(primates),4),]
tree<-Rdnainvar(primates)
}
}
\keyword{phylogenetics}
\keyword{inference}

