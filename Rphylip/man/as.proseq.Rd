\name{as.proseq}
\alias{as.proseq}
\alias{as.phylip.data}
\title{Converts objects to protein sequences or phylip.data object}
\usage{
as.proseq(x, ...)
as.phylip.data(x, ...)
}
\arguments{
	\item{x}{an object containing amino sequences, a set of binary characters, or an arbitrary user-defined multistate character. (Presently only objects of class \code{"phyDat"} are permitted.)}
	\item{...}{optional arguments.}
}
\description{
	Converts objects to class \code{"proseq"} or \code{"phylip.data"}.
}
\value{
	An object of class \code{"proseq"} containing protein sequences; or an object of class \code{"phylip.data"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{print.proseq}}, \code{\link{Rproml}}
}
\keyword{utilities}
\keyword{amino acid}
\keyword{binary character}

