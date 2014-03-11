\name{Rconsense}
\alias{Rconsense}
\title{R interface for consense}
\usage{
Rconsense(trees, path=NULL, ...)
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"}.}
	\item{path}{path to the directory containing the executable consense. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to consense. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{method} which can be \code{"extended"} (extended majority rule consensus, the default), \code{"strict"} (strict consensus), or regular majority rule consensus (\code{"majority"}); \code{"outgroup"} single taxon label or vector of taxa that should be used to root all trees before analysis; \code{rooted} logical value indicated whether to treat the input trees as rooted (defaults to \code{rooted = FALSE}); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the consense program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/consense.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for consense in the PHYLIP package (Felsenstein 2013). consense can be used to compute the consensus tree from a set of phylogenies. 
}
\value{
	This function returns an object of class \code{"phylo"}. For methods other than \code{method = "strict"}, \code{tree$node.label} contains the proportion of phylogenies in \code{trees} containing that subtree.
}
\references{
	Margush, T., McMorris, F.R. (1981) Consensus n-trees. \emph{Bulletin of Mathematical Biology}, 43, 239-244.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\examples{
\dontrun{
trees<-rmtree(n=10,N=10)
tree<-Rconsense(trees)
}
}
\keyword{phylogenetics}
\keyword{consensus}
