source("Rphylip.R")

## sim code for Rcontrast

tree<-pbtree(n=26,scale=1,tip.label=LETTERS)
vcv<-matrix(c(1,1.2,0,1.2,2,0,0,0,3),3,3)
X<-sim.corrs(tree,vcv)
rr<-Rcontrast(tree,X)

n<-round(runif(length(tree$tip.label),1,4))
X<-cbind(sampleFrom(X[,1],xvar=0.1,n=n),sampleFrom(X[,2],xvar=0.1,n=n),sampleFrom(X[,3],xvar=0.1,n=n))
qq<-Rcontrast(tree,X)