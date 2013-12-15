\name{Rdnapenny}
\alias{Rdnapenny}
\title{R interface for dnapenny}
\usage{
Rdnapenny(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnapenny. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to dnapenny. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{groups} number of groups of 1,000 trees (defaults to \code{groups = 10000}); \code{report} reporting frequency, in numbers of trees (defaults to \code{report = 1000}); \code{simple} simple branch & bound (defaults to \code{simple = TRUE}); \code{threshold} threshold value for threshold parsimony (defaults to ordinary parsimony); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the dnapenny program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/dnapenny.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for dnapenny in the PHYLIP package (Felsenstein 2013). dnapenny performs branch & bound parsimony searching following Hendy & Penny (1982).
}
\value{
	This function returns an object of class \code{"phylo"} or \code{"multiPhylo"} that is the tree or trees with the best parsimony score. \code{tree$score} gives the parsimony score, for \code{"phylo"} object \code{tree}.
}
\references{
	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Hendy, M.D., Penny, D. (1982) Branch and bound algorithms to determine minimal evolutionary trees. \emph{Mathematical Biosciences}, 60, 133-142.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnapars}}
}
\examples{
\dontrun{
data(primates)
tree<-Rdnapenny(primates)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
