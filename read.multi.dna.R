read.multi.dna<-function(file,N){
	X<-list()
	skip<-0
	for(i in 1:N){
		if(i==1) X[[i]]<-read.dna(file,format="sequential")
		else {
			skip<-skip+nrow(X[[i-1]])+1
			X[[i]]<-read.dna(file,format="sequential",skip=skip)
		}
	}
	return(X)
}

library(phangorn)
trees<-lapply(X,function(x) NJ(dist.dna(x))) 
class(trees)<-"multiPhylo"