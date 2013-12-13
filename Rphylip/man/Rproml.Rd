\name{Rproml}
\alias{Rproml}
\alias{Rpromlk}
\title{R interfaces for proml and promlk}
\usage{
Rproml(X, path=NULL, ...)
Rpromlk(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"proseq"}.}
	\item{path}{path to the executable containing proml. If \code{path = NULL}, the R will search several commonly used directories for the correct executable file.}
	\item{...}{optional arguments to be passed to proml or promlk. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{tree} object of class \code{"phylo"} - if supplied, then the model will be optimized on a fixed input topology; \code{model} amino acid model - could be \code{"JTT"} (Jones et al. 1992), \code{"PMB"} (Veerassamy et al. 2003), or \code{"PAM"} (Dayhoff & Eck 1968; Dayhoff et al. 1979; Koisol & Goldman 2005); \code{rates} vector of rates (defaults to single rate); \code{rate.categories} vector of rate categories corresponding to the order of \code{rates}; \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity); \code{ncat} number of rate categories for the gamma model; \code{inv} proportion of invariant sites for the invariant sites model (defaults to \code{inv = 0}); \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); \code{speedier} speedier but rougher analysis (defaults to \code{speedier = FALSE}); \code{global} perform global search (defaults to \code{global = TRUE}); \code{random.order} add taxa to tree in random order (defaults to \code{random.order = TRUE}); \code{random.addition} number of random addition replicates for \code{random.order = TRUE} (defaults to \code{random.addition = 10}); \code{outgroup} outgroup if outgroup rooting of the estimated tree is desired; and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	Finally \code{clock=TRUE} enforces a molecular clock. The argument \code{clock} is only available for \code{Rproml}. If \code{clock=TRUE} then promlk is used internally. For \code{Rpromlk} a molecular clock is assumed, thus \code{Rproml(...,clock=TRUE)} and \code{Rpromlk(...)} are equivalent. Note that in PHYLIP 3.695 my tests of promlk yielded peculiar results (all branch lengths zero length, random topology), so I'm not sure what to make of that.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 2013) should first be installed. More information about installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for proml in the PHYLIP package (Felsenstein 1989, 2013). proml can be used for ML phylogeny estimation from amino acid sequences.
}
\value{
	This function returns an object of class \code{"phylo"} that is the optimized tree.
}
\references{
	Dayhoff, M.O., Eck, R.V. (1968) \emph{Atlas of Protein Sequence and Structure 1967-1968}. National Biomedical Research Foundation, Silver Spring, Maryland.

	Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1979) A model of evolutionary change in proteins. pp. 345-352 in \emph{Atlas of Protein Sequence and Structure, Volume 5, Supplement 3, 1978}, ed. M.O. Dayhoff. National Biomedical Research Foundataion, Silver Spring, Maryland.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	Jones, D.T., Taylor, W.R., Thornton, J.M. (1992) The rapid generation of mutation data matrices from protein sequences. \emph{Computer Applications in the Biosciences (CABIOS)}, 8, 275-282.

	Koisol, C., Goldman, N. (2005) Different versions of the Dayhoff rate matrix. \emph{Molecular Biology and Evolution}, 22, 193-199.

	Veerassamy, S., Smith, A., Tillier, E.R. (2003) A transition probability model for amino acid substitutions from blocks. \emph{Journal of Computational Biology}, 10, 997-1010.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{as.proseq}}, \code{\link{Rdnaml}}, \code{\link{read.protein}}
}
\examples{
\dontrun{
data(chloroplast)
tree<-Rproml(chloroplast)
}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
\keyword{amino acid}
