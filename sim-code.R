source("Rphylip.R")

## sim code for Rthreshml

tree<-pbtree(n=10)
X<-fastBM(tree,nsim=4)
Y<-data.frame(X[,1],sapply(X[,2],threshState,th),X[,3],sapply(X[,4],threshState,th))
names(Y)<-paste("X",1:4,sep="")
X<-Y


## sim code for Rcontrast

tree<-pbtree(n=26,scale=1,tip.label=LETTERS)
vcv<-matrix(c(1,1.2,0,1.2,2,0,0,0,3),3,3)
X<-sim.corrs(tree,vcv)
rr<-Rcontrast(tree,X)

n<-round(runif(length(tree$tip.label),1,4))
X<-cbind(sampleFrom(X[,1],xvar=0.1,n=n),sampleFrom(X[,2],xvar=0.1,n=n),sampleFrom(X[,3],xvar=0.1,n=n))
qq<-Rcontrast(tree,X)