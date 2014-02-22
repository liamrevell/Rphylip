\name{print.proseq}
\alias{print.proseq}
\alias{print.phylip.data}
\alias{print.rest.data}
\title{Print method protein sequences, phylip.data, or rest.data object}
\usage{
\method{print}{proseq}(x, printlen=6, digits=3, ...)
}
\arguments{
	\item{x}{an object of class \code{"proseq"}, an object of class \code{"phylip.data"}, or an object of class \code{"rest.data"}.}
	\item{printlen}{number of sequence names to print.}
	\item{digits}{number of digits to print.}
	\item{...}{optional arguments.}
}
\description{
	Print method for an objects of class \code{"proseq"}, \code{"phylip.data"}, or \code{"rest.data"}.
}
\value{
	Prints to screen.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{as.proseq}}, \code{\link{print.DNAbin}}, \code{\link{Rproml}}
}
\keyword{utilities}
\keyword{amino acid}
\keyword{binary character}
\keyword{restriction sites}

