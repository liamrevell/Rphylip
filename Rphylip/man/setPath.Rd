\name{setPath}
\alias{setPath}
\alias{clearPath}
\title{Set path to the folder containing PHYLIP executables for the current R session}
\usage{
setPath(path)
clearPath()
}
\arguments{
	\item{path}{path to the folder contains PHYLIP executables.}
}
\description{
	Sets the path to the folder containing PHYLIP executables for the current R session. Once \code{setPath} has been used, Rphylip functions will no longer search for the path to the PHYLIP executable, using \code{path} instead.

	\code{clearPath} clears the path.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{utilities}

