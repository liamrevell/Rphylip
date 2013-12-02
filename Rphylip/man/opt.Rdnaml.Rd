\name{opt.Rdnaml}
\alias{opt.Rdnaml}
\title{Parameter optimizer for Rdnaml}
\usage{
opt.Rdnaml(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnaml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file. (Currently works only for \code{.Platform$OS.type = "windows"}.)}
	\item{...}{optional arguments. See details for more information.}
}
\details{
	Optional arguments include the following: \code{tree} fixed tree to use in optimization - if not provided, it will be estimated using \code{\link{Rdnaml}} under the default conditions.
}
\description{
	This function is an wrapper for \code{\link{Rdnaml}} that attempts to optimize \code{gamma} (the alpha shape parameter of the gamma model of rate heterogeneity among sites), \code{kappa} (the transition:transversion ratio), and \code{bf} (the base frequencies).
}
\value{
	This function returns a list with the following components: \code{kappa}, \code{gamma}, \code{bf} (see Details), and \code{logLik} (the log-likelihood of the fitted model).
}
\references{
	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{opt.Rdnaml}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
