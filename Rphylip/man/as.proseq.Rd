\name{as.proseq}
\alias{as.proseq}
\alias{as.phylip.data}
\alias{as.rest.data}
\title{Converts objects to protein sequences, phylip.data, or rest.data object}
\usage{
as.proseq(x, ...)
as.phylip.data(x, ...)
as.rest.data(x, ...)
}
\arguments{
	\item{x}{an object containing amino sequences, a set of binary characters, an arbitrary user-defined multistate character, or presences (\code{"+"}) and absences (\code{"-"}) of restriction sites. \code{as.proseq} and \code{as.phylip.data} also converts objects of class \code{"phyDat"} from the phangorn package.}
	\item{...}{optional arguments.}
}
\description{
	Converts objects to class \code{"proseq"}, \code{"phylip.data"}, or \code{"rest.data"}.
}
\value{
	An object of class \code{"proseq"} containing protein sequences; an object of class \code{"phylip.data"}; or an object of class \code{"rest.data"}.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\seealso{
	\code{\link{print.proseq}}, \code{\link{Rproml}}
}
\keyword{utilities}
\keyword{amino acid}
\keyword{binary character}
\keyword{restriction sites}

