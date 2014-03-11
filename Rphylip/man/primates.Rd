\name{primates}
\alias{primates}
\alias{chloroplast}
\alias{primates.bin}
\alias{cotton}
\alias{restriction.data}
\title{Example datasets}
\description{
	\code{primates} is an object of class \code{"DNAbin"} containing nucleotide sequence data of mysterious origin for 12 species of primates. \code{chloroplast} is an object of class\code{"proseq"} containing a chloroplast alignment from the phangorn package (Schliep 2011). \code{primates.bin} is a an object of class \code{"phylip.data"} containing \code{"0"} and \code{"1"} character data for primates. \code{cotton} is a list with two parts: \code{cotton$tree} is an object of class \code{"phylo"}, containing a phylogenetic tree for 37 species of cotton (\emph{Gossypium}); \code{cotton$data} is an object of class \code{"matrix"} containing traits for each species. Most trait values are real, but a few values were not available and so random values were used so there were no missing values in the dataset. \code{restriction.data} is an object of class \code{"rest.data"} containing an example dataset of restriction site data.
}
\usage{
data(primates)
data(chloroplast)
data(primates.bin)
data(cotton)
}
\format{
	The data are stored as an object of class \code{"DNAbin"}, \code{"proseq"}, \code{"phylip.data"}, or \code{"rest.data"}. \code{cotton} is a list containing an object \code{"phylo"} and a matrix of continuous trait data.
}
\source{
	Chamberlain, S.A., & J.A. Rudgers. 2012. How do plants balance multiple mutualists? Correlations among traits for attracting protective bodyguards and pollinators in cotton (Gossypium). Evolutionary Ecology 26:65-77. DOI: 10.1007/s10682-011-9497-3.
}
\author{	
	Liam J. Revell, Scott A. Chamberlain

	Maintainer: Liam J. Revell <liam.revell@umb.edu>
}
\keyword{datasets}
