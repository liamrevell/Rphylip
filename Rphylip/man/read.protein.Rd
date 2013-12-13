\name{read.protein}
\alias{read.protein}
\title{Reads protein sequences from file in multiple formats}
\usage{
read.protein(file, format="fasta", ...)
}
\arguments{
	\item{file}{file name for file containing protein sequences.}
	\item{format}{format of input file. Permitted formats are \code{"fasta"} and \code{"sequential"}. See \code{\link{read.dna}} for more information.}
	\item{...}{optional arguments.}
}
\description{
	Reads protein sequences from a file.
}
\value{
	An object of class \code{"proseq"} containing protein sequences.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{as.proseq}}, \code{\link{print.proseq}}, \code{\link{Rproml}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{amino acid}

