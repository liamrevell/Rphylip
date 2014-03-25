\name{Rseqboot}
\alias{Rseqboot}
\title{R interface for seqboot}
\usage{
Rseqboot(X, path=NULL, ...)
}
\arguments{
	\item{X}{either: (a) an object of class \code{"DNAbin"} containing DNA sequences or an object of class \code{"proseq"} containing amino acid sequences; (b) an object of class \code{"phylip.data"} containing data for a binary or multi-state discrete character; (c) an object of class \code{"rest.data"} containing restriction site data; or (d) a matrix with gene frequency data (this option does not presently work properly).}
	\item{path}{path to the executable containing seqboot. If \code{path = NULL}, the R will search several paths for seqboot. See details for more information.}
	\item{...}{optional arguments to be passed to seqboot. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{type} the data type (options are \code{"sequence"}, \code{"morph"}, \code{"rest"}, or \code{"gene.freq"}); \code{method} analysis method (options are \code{"bootstrap"}, \code{"jacknife"}, and \code{"permute"}); \code{percentage} a bootstrap or jacknife percentage (can be \code{>100}); \code{block.size} block size in number of characters for the bootstrap, jacknife, or permutation (defaults to \code{1}, the regular bootstrap); \code{replicates} number of replicates (defaults to \code{100}); \code{weights} a vector of weights that must be \code{0} or \code{1} - characters with weight of \code{0} are excluded from the resampling procedure; \code{rate.categories} vector of integers between \code{1} and \code{9} indicating the rate category for each character - this is supplied so that resampled rate categories can be output by the method; \code{mixture} vector of \code{"W"} and \code{"C"} (see \code{\link{Rmix}}); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the seqboot program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/seqboot.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for seqboot in the PHYLIP package (Felsenstein 2013). seqboot can be to perform the non-parametric bootstrap following Felsenstein (1985).
}
\value{
	This function returns a list of objects of class \code{"DNAbin"}, \code{"proseq"}, \code{"phylip.data"}, or \code{"rest.data"} depending on the input data type. If \code{ancestors}, \code{mixture}, or \code{rate.categories} are supplied, then each element of the list contains an object \code{data} (containing the data), and vectors \code{ancestors}, \code{mixture}, and/or \code{categories}, containing the resampled input vectors.
}
\references{
	Felsenstein, J. (1985) Confidence limits on phylogenies: an approach using the bootstrap. \emph{Evolution} 39, 783-791.  

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle. 
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\keyword{phylogenetics}
