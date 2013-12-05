require(phytools)

source("Rphylip.R")

## sim code for Rtreedist

trees<-pbtree(n=26,tip.label=LETTERS,nsim=7)
trees2<-pbtree(n=26,tip.label=LETTERS,nsim=17)
D<-Rtreedist(trees,trees2=trees2,method="symmetric",quiet=TRUE)
D<-Rtreedist(trees,method="symmetric")
D<-Rtreedist(trees,method="branch.score")


system.time(Rtreedist(trees,method="symmetric"))
system.time(multiRF(trees))

## sim code for Rthreshml

tree<-pbtree(n=100)
x<-fastBM(tree,nsim=2)
th<-setNames(c(0,Inf),c(0,1))
X<-data.frame(x[,1],sapply(x[,1],threshState,th),x[,2])
names(X)<-paste("v",1:2,sep="")
fit<-Rthreshml(tree,X,proposal=0.2)


## sim code for Rcontrast

tree<-pbtree(n=26,scale=1,tip.label=LETTERS)
vcv<-matrix(c(1,1.2,0,1.2,2,0,0,0,3),3,3)
X<-sim.corrs(tree,vcv)
rr<-Rcontrast(tree,X)

n<-round(runif(length(tree$tip.label),1,4))
X<-cbind(sampleFrom(X[,1],xvar=0.1,n=n),sampleFrom(X[,2],xvar=0.1,n=n),sampleFrom(X[,3],xvar=0.1,n=n))
qq<-Rcontrast(tree,X)