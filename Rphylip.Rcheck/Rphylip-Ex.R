pkgname <- "Rphylip"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('Rphylip')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rconsense")
### * Rconsense

flush(stderr()); flush(stdout())

### Name: Rconsense
### Title: R interface for consense
### Aliases: Rconsense
### Keywords: phylogenetics

### ** Examples

## Not run: 
##D trees<-rmtree(n=10,N=10)
##D tree<-Rconsense(trees)
## End(Not run)



cleanEx()
nameEx("Rdnadist")
### * Rdnadist

flush(stderr()); flush(stdout())

### Name: Rdnadist
### Title: R interfaces for dnadist
### Aliases: Rdnadist
### Keywords: phylogenetics inference maximum likelihood distance method

### ** Examples

## Not run: 
##D data(primates)
##D D<-Rdnadist(primates,kappa=10)
##D tree<-Rneighbor(D)
## End(Not run)



cleanEx()
nameEx("Rdnaml")
### * Rdnaml

flush(stderr()); flush(stdout())

### Name: Rdnaml
### Title: R interfaces for dnaml and dnamlk
### Aliases: Rdnaml Rdnamlk
### Keywords: phylogenetics inference maximum likelihood

### ** Examples

## Not run: 
##D data(primates)
##D tree<-Rdnaml(primates,kappa=10)
##D clockTree<-Rdnamlk(primates,kappa=10)
## End(Not run)



cleanEx()
nameEx("Rdnapars")
### * Rdnapars

flush(stderr()); flush(stdout())

### Name: Rdnapars
### Title: R interface for dnapars
### Aliases: Rdnapars
### Keywords: phylogenetics inference parsimony

### ** Examples

## Not run: 
##D data(primates)
##D tree<-Rdnapars(primates)
## End(Not run)



cleanEx()
nameEx("Rdnapenny")
### * Rdnapenny

flush(stderr()); flush(stdout())

### Name: Rdnapenny
### Title: R interface for dnapenny
### Aliases: Rdnapenny
### Keywords: phylogenetics inference parsimony

### ** Examples

## Not run: 
##D data(primates)
##D tree<-Rdnapenny(primates)
## End(Not run)



cleanEx()
nameEx("Rneighbor")
### * Rneighbor

flush(stderr()); flush(stdout())

### Name: Rneighbor
### Title: R interface for neighbor
### Aliases: Rneighbor
### Keywords: phylogenetics inference distance method

### ** Examples

## Not run: 
##D data(primates)
##D D<-dist.dna(data(primates),model="JC")
##D tree<-Rneighbor(D)
## End(Not run)



cleanEx()
nameEx("Rproml")
### * Rproml

flush(stderr()); flush(stdout())

### Name: Rproml
### Title: R interfaces for proml and promlk
### Aliases: Rproml Rpromlk
### Keywords: phylogenetics inference maximum likelihood

### ** Examples

## Not run: 
##D data(chloroplast)
##D tree<-Rproml(chloroplast)
## End(Not run)



cleanEx()
nameEx("Rtreedist")
### * Rtreedist

flush(stderr()); flush(stdout())

### Name: Rtreedist
### Title: R interface for treedist
### Aliases: Rtreedist
### Keywords: phylogenetics

### ** Examples

## Not run: 
##D trees<-rmtree(n=10,N=10)
##D D<-Rtreedist(trees,method="symmetric")
## End(Not run)



cleanEx()
nameEx("opt.Rdnaml")
### * opt.Rdnaml

flush(stderr()); flush(stdout())

### Name: opt.Rdnaml
### Title: Parameter optimizer for Rdnaml
### Aliases: opt.Rdnaml
### Keywords: phylogenetics inference maximum likelihood

### ** Examples

## Not run: 
##D data(primates)
##D fit<-opt.Rdnaml(primates,bounds=list(kappa=c(0.1,40))
##D tree<-Rdnaml(primates,kappa=fit$kappa,gamma=fit$gamma,bf=fit$bf)
## End(Not run)



cleanEx()
nameEx("setupOSX")
### * setupOSX

flush(stderr()); flush(stdout())

### Name: setupOSX
### Title: Help set up PHYLIP in Mac OS X
### Aliases: setupOSX
### Keywords: phylogenetics utilities

### ** Examples

## Not run: 
##D setupOSX()
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
