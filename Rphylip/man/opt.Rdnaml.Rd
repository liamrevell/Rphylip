\name{opt.Rdnaml}
\alias{opt.Rdnaml}
\title{Parameter optimizer for Rdnaml}
\usage{
opt.Rdnaml(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"DNAbin"}.}
	\item{path}{path to the executable containing dnaml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments. See details for more information.}
}
\details{
	Optional arguments include the following: \code{tree} fixed tree to use in optimization - if not provided, it will be estimated using \code{\link{Rdnaml}} under the default conditions; \code{bounds} a list with bounds for optimization - for \code{kappa} and \code{gamma} this should be a two-element vector, whereas for \code{bf} this should be a 4 x 2 matrix with lower bounds in column 1 and upper bounds in column 2.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an wrapper for \code{\link{Rdnaml}} that attempts to optimize \code{gamma} (the alpha shape parameter of the gamma model of rate heterogeneity among sites), \code{kappa} (the transition:transversion ratio), and \code{bf} (the base frequencies).
}
\value{
	This function returns a list with the following components: \code{kappa}, \code{gamma}, \code{bf} (see Details), and \code{logLik} (the log-likelihood of the fitted model).
}
\references{
	Felsenstein, J. (1981) Evolutionary trees from DNA sequences: A Maximum Likelihood approach. \emph{Journal of Molecular Evolution}, 17, 368-376.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Felsenstein, J., Churchill, G. A. (1996) A Hidden Markov Model approach to variation among sites in rate of evolution. \emph{Molecular Biology and Evolution}, 13, 93-104.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Rdnaml}}
}
\examples{
\dontrun{
data(primates)
fit<-opt.Rdnaml(primates,bounds=list(kappa=c(0.1,40))
tree<-Rdnaml(primates,kappa=fit$kappa,gamma=fit$gamma,bf=fit$bf)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
\keyword{DNA}
