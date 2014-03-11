\name{Rprotpars}
\alias{Rprotpars}
\title{R interface for protpars}
\usage{
Rprotpars(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"proseq"} containing aligned amino acid sequences.}
	\item{path}{path to the executable containing protpars. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to protpars. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the parsimony score will be computed on a fixed input topology; \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{threshold} threshold value for threshold parsimony (defaults to ordinary parsimony); \code{genetic.code}, type of genetic code to assume (options are \code{"universal"}, the default, \code{"mitochondrial"}, \code{"vertebrate.mitochondrial"}, \code{"fly.mitochondrial"}, and \code{"yeast.mitochondrial"}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the protpars program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/protpars.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for protpars in the PHYLIP package (Felsenstein 2013). protpars can be used for MP phylogeny estimation from protein sequences (Eck & Dayhoff 1966; Fitch 1971).
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Eck. R.V., Dayhoff, M.O. (1966) \emph{Atlas of Protein Sequence and Structure 1966}. National Biomedical Research Foundation, Silver Spring, Maryland.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Fitch, W.M. (1971) Toward defining the course of evolution: Minimu change for a specified tree topology. \emph{Systematic Zoology}, 20, 406-416.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\seealso{
	\code{\link{as.proseq}}, \code{\link{Rdnapars}}, \code{\link{read.protein}}
}
\examples{
\dontrun{
data(chloroplast)
tree<-Rprotpars(chloroplast)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
\keyword{amino acid}
