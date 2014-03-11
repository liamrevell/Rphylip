\name{Rmix}
\alias{Rmix}
\title{R interface for mix}
\usage{
Rmix(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"phylip.data"} with a binary character; or a matrix. Data should be \code{"0"}, \code{"1"}, and \code{"?"}.}
	\item{path}{path to the executable containing mix. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to mix. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the parsimony score will be computed on a fixed input topology; \code{mixture}, a vector containing \code{"W"} and \code{"C"} in upper or lower case indicating which characters should be evaluated using Wagner (Eck & Dayhoff 1966; Kluge & Farris 1969) or Camin-Sokal (Camin & Sokal 1965) parsimony methods, respectively - this argument takes precedence over \code{method} (see below); \code{method} either \code{"Wagner"} or \code{"Camin-Sokal"}, in upper or lower case, indicating whether Wagner (Eck & Dayhoff 1966; Kluge & Farris 1969) or Camin-Sokal (Camin & Sokal 1965) parsimony should be used; \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{threshold} threshold value for threshold parsimony (defaults to ordinary parsimony); \code{ancestral} a vector of length equal to the number of columns in \code{X} containing the the ancestral state at the root of the tree for each character, these should be \code{"0"}, \code{"1"}, and \code{"?"}; \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the mix program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/mix.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for mix in the PHYLIP package (Felsenstein 2013). mix performs branch & bound parsimony searching following Hendy & Penny (1982).
}
\value{
	This function returns an object of class \code{"phylo"} or \code{"multiPhylo"} that is the tree or trees with the best parsimony score. \code{tree$score} gives the parsimony score, for \code{"phylo"} object \code{tree}.
}
\references{
	Camin, J.H., Sokal, R.R. (1965) A method for deducing branching sequences in phylogeny. \emph{Evolution}, 19, 311-326.

	Eck, R.V., Dayhoff, M.O. (1966) \emph{Atlas of Protein Sequence and Structure 1966}. National Biomedical Research Foundation, Silver Spring, Maryland.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Kluge, A.G., Farris, J.S. (1969) Quantitative phyletics and the evolution of anurans. \emph{Systematic Zoology}, 18, 1-32.

	Hendy, M.D., Penny, D. (1982) Branch and bound algorithms to determine minimal evolutionary trees. \emph{Mathematical Biosciences}, 60, 133-142.
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
data(primates.bin)
tree<-Rmix(primates.bin)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
