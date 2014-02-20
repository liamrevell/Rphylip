\name{Rfitch}
\alias{Rfitch}
\title{R interface for fitch}
\usage{
Rfitch(D, path=NULL , ...)
}
\arguments{
	\item{D}{a distance matrix as an object of class \code{"matrix"} or \code{"dist"}. If a matrix, then \code{D} should be symmetrical with a diagonal of zeros.}
	\item{path}{path to the executable containing fitch. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to fitch. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{method} - can be \code{"FM"} or \code{"fm"} (for the Fitch-Margoliash criterion; Fitch & Margoliash 1967), \code{"LS"} or \code{"ls"} for the Cavalli-Sforza & Edwards (1967) least-squares criterion, or \code{"ME"} or \code{"me"} (for Minimum Evolution); \code{tree} object of class \code{"phylo"} - if supplied, then branch lengths are optimized & the optimality criterion is computed on the input tree; \code{negative} a logical value indicating whether negative branch lengths should be permitted (defaults to \code{negative = TRUE}; \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the fitch program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/fitch.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for fitch in the PHYLIP package (Felsenstein 2013). fitch can be used for tree inference using the Fitch-Margoliash method (Fitch & Margoliash 1967), the Cavalli-Sforza & Edwards leasts-squares method (Cavalli-Sforza & Edwards 1967), and the minimum evolution criterion (Kidd & Sgaramella-Zonta 1971; Nei & Rzhetsky 1993).
}
\value{
	This function returns an object of class \code{"phylo"} that is the Fitch-Margoliash, LS, or ME tree.
}
\references{
	 Cavalli-Sforza, L.L., Edwards, A.W.F. (1967) Phylogenetic analysis: models and estimation procedures. \emph{Evolution}, 32, 550-570. 

	 Fitch, W.M., Margoliash, E. (1967) Construction of phylogenetic trees. \emph{Science}, 155, 279-284.

	 Kidd, K.K., Sgaramella-Zonta, L.A. (1971) Phylogenetic analysis: concepts and methods. \emph{American Journal of Human Genetics}, 23, 235-252. 	

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnadist}}
}
\examples{
\dontrun{
data(primates)
D<-dist.dna(data(primates),model="JC")
tree<-Rfitch(D)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{distance method}
