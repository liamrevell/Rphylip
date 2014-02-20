\name{Rgendist}
\alias{Rgendist}
\title{R interface for gendist}
\usage{
Rgendist(X, path=NULL, ...)
}
\arguments{
	\item{X}{either: (a) A \emph{matrix} of allele frequencies (in columns) with rownames containing species names. In this case it is assumed that \code{m - 1} allele frequencies are reported for a locus with \code{m} segregating alleles. (The \code{m}th frequency is just one minus the sum of the other \code{m-1}.) Or (b) a list of matrices in which each row contains the relative frequency of alleles at a locus for a species. In this case, the rownames of each matrix in the list should contain the species names. Here the number of alleles at each locus can be inferred from the number of columns in the matrix, and no frequencies should be omitted. (I.e., rows should add to 1.0.)}
	\item{path}{path to the executable containing gendist. If \code{path = NULL}, the R will search several paths for gendist. See details for more information.}
	\item{...}{optional arguments to be passed to gendist. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{method}, the method for calculating the genetic distance - options are \code{"Nei"}, \code{"Cavalli-Sforza"}, and \code{"Reynolds"} (defaults to \code{method="Nei"}); \code{nalleles} a vector containing the number of alleles per locus in \code{X} (not used if \code{X} is a list, and assumed to be a vector of \code{2}s if ommitted); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the gendist program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/gendist.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for gendist in the PHYLIP package (Felsenstein 2013). gendist can be used to estimate the evolutionary distances between populations based on gene frequency data. Three methods can be used for computing these distances. These are (1) Nei's genetic distance (Nei 1972), (2) Cavalli-Sforza's chord measure (Cavalli-Sforza & Edwards 1967), and (3) Reynolds et al.'s (1983) genetic distance.
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{
	Cavalli-Sforza, L.L., Edwards, A.W.F. (1967) Phylogenetic analysis: models and estimation procedures. \emph{Evolution}, 32, 550-570. 

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Nei, M. (1972) Genetic distance between populations. \emph{American Naturalist}, 106, 283-292. 

	Reynolds, J.B., Weir, B.S., Cockerham, C.C. (1983) Estimation of the coancestry coefficient: basis for a short-term genetic distance. \emph{Genetics}, 105, 767-779. 
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnadist}}, \code{\link{Rprotdist}}
}
\examples{

}
\keyword{phylogenetics}
\keyword{distance method}
