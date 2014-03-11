\name{Rprotdist}
\alias{Rprotdist}
\title{R interfaces for protdist}
\usage{
Rprotdist(X, path=NULL, ...)
}
\arguments{
	\item{X}{an object of class \code{"proseq"} containing aligned amino acid sequences.}
	\item{path}{path to the executable containing protdist. If \code{path = NULL}, the R will search sev arguments to be passed to protdist. See details for more information.}
	\item{...}{optional arguments to be passed to protdist. See details for more information.}
}
\details{
	Optional arguments include the following: \code{quiet} suppress some output to R console (defaults to \code{quiet = FALSE}); \code{model} can be \code{"JTT"} (Jones et al. 1992), \code{"PMB"} (Veerassamy et al. 2003), \code{"PAM"} (Dayhoff & Eck 1968; Dayhoff et al. 1979; Koisol & Goldman 2005), \code{"Kimura"} (a simple model based on Kimura 1980), \code{"similarity"} which gives the similarity between sequences, and \code{"categories"} which is due to Felsenstein; \code{gamma} alpha shape parameter of a gamma model of rate heterogeneity among sites (defaults to no gamma rate heterogeneity) - note that gamma rate heterogeneity does not apply to \code{model = "Kimura"} or \code{model = "similarity"}; \code{kappa} transition:transversion ratio (defaults to \code{kappa = 2.0}), \code{genetic.code}, type of genetic code to assume (options are \code{"universal"}, the default, \code{"mitochondrial"}, \code{"vertebrate.mitochondrial"}, \code{"fly.mitochondrial"}, and \code{"yeast.mitochondrial"}), \code{categorization}, categorization scheme for amino acids (options are \code{"GHB"}, the George et al. 1988 classification, \code{"Hall"}, a classification scheme provided by Ben Hall, and \code{"chemical"}, a scheme based on Conn & Stumpf 1963); and, finally, \code{ease}, a numerical parameter that indicates the facility of getting between amino acids of different categories in which 0 is nearly impossible, and 1 is no difficulty (defaults to \code{ease = 0.457}) - note that \code{kappa}, \code{bf}, \code{genetic.code}, \code{categorization}, and \code{ease} apply only to \code{model = "categories"}; \code{rates} vector of rates (defaults to single rate); \code{rate.categories} vector of rate categories corresponding to the order of \code{rates}; \code{weights} vector of weights of length equal to the number of columns in \code{X} (defaults to unweighted); and \code{cleanup} remove PHYLIP input & output files after the analysis is completed (defaults to \code{cleanup = TRUE}).

	More information about the protdist program in PHYLIP can be found here \url{http://evolution.genetics.washington.edu/phylip/doc/protdist.html}.

	Obviously, use of any of the functions of this package requires that PHYLIP (Felsenstein 1989, 2013) should first be installed. Instructions for installing PHYLIP can be found on the PHYLIP webpage: \url{http://evolution.genetics.washington.edu/phylip.html}.
}
\description{
	This function is an R interface for protdist in the PHYLIP package (Felsenstein 2013). protdist can be used to estimate the evolutionary distances between amino acid sequences under various models.
}
\value{
	This function returns an object of class \code{"dist"}.
}
\references{
	 Conn, E.E., Stumpf, P.K. (1963) \emph{Outlines of Biochemistry}. John Wiley and Sons, New York. 

	Dayhoff, M.O., Eck, R.V. (1968) \emph{Atlas of Protein Sequence and Structure 1967-1968}. National Biomedical Research Foundation, Silver Spring, Maryland.

	Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1979) A model of evolutionary change in proteins. pp. 345-352 in \emph{Atlas of Protein Sequence and Structure, Volume 5, Supplement 3, 1978}, ed. M.O. Dayhoff. National Biomedical Research Foundataion, Silver Spring, Maryland.

	Felsenstein, J. (1989) PHYLIP--Phylogeny Inference Package (Version 3.2). \emph{Cladistics}, 5, 164-166.

	Felsenstein, J. (2013) PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

	 George, D.G., Hunt, L.T., Barker., W.C. (1988) \emph{Current methods in sequence comparison and analysis}. pp. 127-149 in Macromolecular Sequencing and Synthesis, ed. D. H. Schlesinger. Alan R. Liss, New York. 

	Jones, D.T., Taylor, W.R., Thornton, J.M. (1992) The rapid generation of mutation data matrices from protein sequences. \emph{Computer Applications in the Biosciences (CABIOS)}, 8, 275-282.

	Kimura, M. (1980) A simple model for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. \emph{Journal of Molecular Evolution}, 16, 111-120.

	Koisol, C., Goldman, N. (2005) Different versions of the Dayhoff rate matrix. \emph{Molecular Biology and Evolution}, 22, 193-199.

	Veerassamy, S., Smith, A., Tillier, E.R. (2003) A transition probability model for amino acid substitutions from blocks. \emph{Journal of Computational Biology}, 10, 997-1010.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\seealso{
	\code{\link{Rneighbor}}
}
\examples{
\dontrun{
data(chloroplast)
D<-Rprotdist(chloroplast,model="PAM")
tree<-Rneighbor(D)
}
}
\keyword{phylogenetics}
\keyword{maximum likelihood}
\keyword{distance method}
\keyword{amino acid}
