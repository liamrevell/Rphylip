pkgname <- "Rphylip"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('Rphylip')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rdnaml")
### * Rdnaml

flush(stderr()); flush(stdout())

### Name: Rdnaml
### Title: R interface for dnaml
### Aliases: Rdnaml
### Keywords: phylogenetics inference maximum likelihood

### ** Examples

## Not run: 
##D data(primates)
##D tree<-Rdnaml(primates,kappa=10)
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