\name{setupOSX}
\alias{setupOSX}
\title{Help set up PHYLIP in Mac OS X}
\usage{
setupOSX(path=NULL)
}
\arguments{
	\item{path}{path to the folder containing the PHYLIP package. If \code{path = NULL}, the R will search several commonly used directories.}
}
\details{
	This function can be used to help set up PHYLIP (\url{http://evolution.genetics.washington.edu/phylip.html}) following the special instructions found here: \url{http://evolution.genetics.washington.edu/phylip/install.html}. \code{setupOSX} should only be run once - when PHYLIP is first installed.
}
\description{
	This function attempts to help set up PHYLIP on a Mac OS X machine.
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\examples{
\dontrun{
setupOSX()
}
}
\keyword{phylogenetics}
\keyword{utilities}
