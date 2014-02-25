\name{Rrestml}
\alias{Rrestml}
\title{R interface for restml}
\usage{
Rrestml(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"rest.data"} containing data for restriction sites.}
	\item{path}{path to the executable containing restml. If \code{path = NULL}, the R will search several paths for restml. See details for more information.}
	\item{...}{optional arguments to be passed to restml. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the model will be optimized on a fixed input topology; \code{speedier} speedier but rougher analysis (defaults to \code{speedier = FALSE}); \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{site.length} the restriction site length (defaults to \code{site.length=6}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the restml program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/restml.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for restml in the PHYLIP package (Felsenstein 2013). restml can be used to conduct ML phylogeny inference from restriction site or fragment data (Nei & Li 1979; Smouse & Li 1987; Felsenstein 1992).
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{
	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (1992) Phylogenies from restriction sites, a maximum likelihood approach. \emph{Evolution}, 46, 159-173. 

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Nei, M., Li, W.-H. (1979) Mathematical model for studying genetic variation in terms of restriction endonucleases. \emph{Proceedings of the National Academy of Sciences}, 76, 5269-5273.

	Smouse, P. E., Li, W.-H. (1987) Likelihood analysis of mitochondrial restriction-cleavage patterns for the human-chimpanzee-gorilla trichotomy. \emph{Evolution}, 41, 1162-1176.  
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rrestdist}}
}
\examples{
\dontrun{
data(restriction.data)
tree<-Rrestml(restriction.data)
}
}
\keyword{phylogenetics}
\keyword{distance method}
