## function check before overwriting files
## written by Liam J. Revell 2013

file.warn<-function(gg){
	ff<-list.files()
	gg[sapply(gg,"%in%",ff)]->gg
	if(any(sapply(gg,"%in%",ff))){
		cat(paste("Warning:\n  One or more of",paste("\"",gg,"\"",sep="",collapse=", "),
		"\n  was found in your current working directory and may be overwritten\n"))
		cat("\nPress ENTER to continue or q to quit: ")
		q<-readLines(n=1)
		if(q=="q"||q=="Q") return(0) else return(1)
	} else return(1)
}

## calls dnapars from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnapars<-function(X,path=".",...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(thorough)) thorough<-list(...)$thorough
	else thorough<-TRUE
	if(!thorough) oo<-c(oo,"s","n")
	if(hasArg(nsave)) nsave<-list(...)$nsave
	else nsave<-10000
	if(nsave!=10000) oo<-c(oo,"v",nsave)
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(transversion)) transversion<-list(...)$transversion
	else transversion<-FALSE
	if(transversion) oo<-c(oo,"n")
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[i,],collapse=""),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnapars",sep=""),input=oo)
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		tree[[i]]$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		if(class(tree)=="phylo") tree<-root(tree,outgroup)
		else if(class(tree)=="multiPhylo"){
			tree<-lapply(tree,root,outgroup=outgroup)
			class(tree)<-"multiPhylo"
		}
		if(!quiet){
			cat("Rooted with the outgroup\n")
			cat("------------------------\n")
			cat(paste(paste(outgroup,collapse=", "),"\n\n"))
		}
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		system("rm infile",show.output.on.console=FALSE)
		system("rm outtree",show.output.on.console=FALSE)
		system("rm outfile",show.output.on.console=FALSE)
		if(!is.null(weights)) system("rm weights",show.output.on.console=FALSE)
		if(intree) system("rm intree",show.output.on.console=FALSE)	
	}
	return(tree)
}

## calls contml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontml<-function(X,path=".",...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(is.matrix(X)){
		## assumes X is a matrix of continuous character data
		write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
		for(i in 1:nrow(X)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
		oo<-c(oo,"c")
		if(hasArg(tree)){
			oo<-c(oo,"u")
			tree<-list(...)$tree
			tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
			write.tree(tree,"intree")
			intree<-TRUE
		} else intree<-FALSE
		if(hasArg(global)) global<-list(...)$global
		else global<-TRUE
		if(global) oo<-c(oo,"g")
		if(hasArg(random.order)) random.order<-list(...)$random.order
		else random.order<-TRUE
		if(random.order){
			if(hasArg(random.addition)) random.addition<-list(...)$random.addition
			else random.addition<-10
			oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
		}
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contml",sep=""),input=oo)
		tree<-read.tree("outtree")
		temp<-readLines("outfile")
		logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
		temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		if(!quiet){
			cat("Translation table\n")
			cat("-----------------\n")
			temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",
				sep="")),y=rownames(X))
			cat("\n")
		}
		tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	} else if(is.list(X)){
		## assumes X is a list of matrices containing gene frequency data
		tips<-rownames(X[[1]])
		X<-lapply(X,function(x,tips) x[tips,],tips=tips)
		write(paste("    ",nrow(X[[1]]),"   ",length(X),sep=""),file="infile")
		nalleles<-sapply(X,ncol)
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		## verify that all rows of all X sum to 1.0
		temp<-sapply(X,rowSums)
		if(!all(round(temp,2)==1)) stop("Some of the rows of X do not sum to 1.0")
		for(i in 1:length(tips)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			dd<-vector()
			for(j in 1:length(X)) dd<-c(dd,X[[j]][i,])
			tt<-paste(sp,paste(dd,collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
		oo<-c(oo,"a")
		if(hasArg(tree)){
			oo<-c(oo,"u")
			tree<-list(...)$tree
			tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
			write.tree(tree,"intree")
			intree<-TRUE
		} else intree<-FALSE
		if(hasArg(global)) global<-list(...)$global
		else global<-TRUE
		if(global) oo<-c(oo,"g")
		if(hasArg(random.order)) random.order<-list(...)$random.order
		else random.order<-TRUE
		if(random.order){
			if(hasArg(random.addition)) random.addition<-list(...)$random.addition
			else random.addition<-10
			oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
		}
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contml",sep=""),input=oo)
		tree<-read.tree("outtree")
		temp<-readLines("outfile")
		logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
		temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		if(!quiet){
			cat("Translation table\n")
			cat("-----------------\n")
			temp<-lapply(1:length(tips),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",
				sep="")),y=tips)
			cat("\n")
		}
		tree$tip.label<-tips[as.numeric(tree$tip.label)]
	} else stop("X should be a matrix (for continuous characters) or a list (for gene frequencies)")
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-root(tree,outgroup)
		if(!quiet){
			cat("Rooted with the outgroup\n")
			cat("------------------------\n")
			cat(paste(paste(outgroup,collapse=", "),"\n\n"))
		}
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		system("rm infile",show.output.on.console=FALSE)
		system("rm outtree",show.output.on.console=FALSE)
		system("rm outfile",show.output.on.console=FALSE)
		if(intree) system("rm intree",show.output.on.console=FALSE)	
	}
	tree$logLik<-logLik
	return(tree)
}

## calls dnaml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnaml<-function(X,path=".",...){
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("categories","infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(kappa)){
		kappa<-list(...)$kappa
		oo<-c(oo,"t",kappa)
	}
	if(hasArg(bf)){
		bf<-list(...)$bf
		bf<-bf/sum(bf)
		bf<-paste(bf,collapse=" ")
		oo<-c(oo,"f",bf)
	}
	if(hasArg(rates)){
		rates<-list(...)$rates
		if(hasArg(rate.categories)){
			rate.categories<-list(...)$rate.categories
			write(paste(rate.categories,collapse=""),file="categories")
			ncats<-length(rates)
			rates<-paste(rates,collapse=" ")
			oo<-c(oo,"c",ncats,rates)
		} else {
			warning("cannot use rates argument without rate categories; ignoring argument rates")
			rates<-NULL
		}
	} else rates<-NULL
	if(hasArg(gamma)) gamma<-list(...)$gamma
	else gamma<-NULL
	if(hasArg(inv)) inv<-list(...)$inv
	else inv<-NULL
	if(hasArg(ncat)) ncat<-list(...)$ncat
	else ncat<-4
	if(!is.null(gamma)&&is.null(inv)){
		oo<-c(oo,"r")
		ee<-c(ee,1/sqrt(gamma),ncat)
	} else if(!is.null(gamma)&&!is.null(inv)){
		oo<-c(oo,"r","r")
		ee<-c(ee,1/sqrt(gamma),inv)
	}
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(speedier)) speedier<-list(...)$speedier
	else speedier<-FALSE
	if(!speedier) oo<-c(oo,"s")
	if(hasArg(global)) global<-list(...)$global
	else global<-TRUE
	if(global) oo<-c(oo,"g")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y",ee,"r")
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[i,],collapse=""),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
	system("touch outtree")
	system("touch outfile")
	temp<-system(paste(path,"/dnaml",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-root(tree,outgroup)
		if(!quiet){
			cat("Rooted with the outgroup\n")
			cat("------------------------\n")
			cat(paste(paste(outgroup,collapse=", "),"\n\n"))
		}
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		system("rm infile",show.output.on.console=FALSE)
		system("rm outtree",show.output.on.console=FALSE)
		system("rm outfile",show.output.on.console=FALSE)
		if(!is.null(weights)) system("rm weights",show.output.on.console=FALSE)
		if(!is.null(rates)) system("rm categories",show.output.on.console=FALSE)
		if(intree) system("rm intree",show.output.on.console=FALSE)	
	}
	tree$logLik<-logLik
	return(tree)
}

## function to optimize parameters of Rdnaml
## written by Liam J. Revell 2013

opt.Rdnaml<-function(X,path=".",...){
	if(hasArg(tree)) tree<-list(...)$tree
	else {
		cat("\nFinding starting tree for parameter optimization\n")
		cat("------------------------------------------------\n")
		tree<-Rdnaml(X,path)
	}
	lik<-function(par,X,tree,path){
		kappa<-par[1]
		gamma<-par[2]
		bf<-par[c(3,4,5,6)]/sum(par[c(3,4,5,6)])
		ll<-Rdnaml(X,path,tree=tree,kappa=kappa,gamma=gamma,bf=bf,speedier=TRUE,random.order=FALSE,
			global=FALSE,quiet=TRUE)$logLik
		cat("\nOptimization progress\n")
		cat("-----------------------\n")
		cat(paste("kappa: ",round(kappa,5),"; gamma: ",round(gamma,5),"; bf: [",
			paste(round(bf,4),collapse=","),"]; logLik: ",round(ll,4),"\n",sep=""))
		cat("-----------------------\n")
		ll
	}
	par<-c(1,10,rep(0.25,4))
	fit<-optim(par,lik,X=X,tree=tree,path=path,method="L-BFGS-B",
		lower=c(rep(0.01,2),rep(0.01,4)),upper=c(20,20,rep(1,4)),
		control=list(fnscale=-1))
	return(list(kappa=fit$par[1],gamma=fit$par[2],
		bf=fit$par[3:6]/sum(fit$par[3:6]),logLik=fit$value))
}

## function calls contrast from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontrast<-function(tree,X,path=".",...){
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	if(!is.binary.tree(tree)){
		cat("Warning:\n  Tree is not binary, resolving with branches of zero length\n")
		tree<-multi2di(tree)
	}
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile"))==0) return(NULL)
	oo<-c("r")
	if(length(unique(rownames(X)))==nrow(X)){
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y)[1],y=rownames(X))
		write.tree(tree,"intree")
		write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
		for(i in 1:nrow(X)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
		oo<-c(oo,"c")
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contrast",sep=""),input=oo)
		temp<-readLines("outfile")
		ii<-grep("Contrasts",temp)
		Contrasts<-matrix(NA,tree$Nnode,ncol(X))
		for(i in 1:tree$Nnode){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Contrasts[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Covariance",temp)
		Covariance_matrix<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Covariance_matrix[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Regressions",temp)
		Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Correlations",temp)
		Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Correlations[i,]<-as.numeric(x[x!=""])
		}
		if(hasArg(cleanup)) cleanup<-list(...)$cleanup
		else cleanup<-TRUE
		if(cleanup){
			system("rm infile",show.output.on.console=FALSE)
			system("rm intree",show.output.on.console=FALSE)
			system("rm outfile",show.output.on.console=FALSE)
		}
		if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		return(list(Contrasts=Contrasts,Covariance_matrix=Covariance_matrix,
			Regressions=Regressions,Correlations=Correlations))	
	} else {
		tips<-tree$tip.label
		tree$tip.label<-1:length(tree$tip.label)
		write.tree(tree,"intree")
		write(paste("    ",length(tips),"   ",ncol(X),sep=""),file="infile")
		for(i in 1:length(tree$tip.label)){
			ii<-which(rownames(X)==tips[i])
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,length(ii),collapse=" ")
			write(tt,append=TRUE,file="infile")
			for(j in 1:length(ii)) write(paste(X[ii[j],],collapse=" "),append=TRUE,file="infile")
		}
		oo<-c(oo,"w")
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contrast",sep=""),input=oo)
		temp<-readLines("outfile")
		ii<-grep("Estimate of VarA",temp)
		VarA<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Estimate of VarE",temp)[1]
		VarE<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarA Regressions",temp)
		VarA.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarA Correlations",temp)
		VarA.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[1]
		VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[1]
		VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Correlations",temp)[1]
		VarE.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Estimate of VarE",temp)[2]
		nonVa.VarE<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[2]
		nonVa.VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Correlations",temp)[2]
		nonVa.VarE.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Log likelihood with varA",temp)
		aa<-strsplit(strsplit(temp[ii],"=")[[1]][2]," ")[[1]]
		aa<-aa[aa!=""]
		logLik<-as.numeric(sub(",","",aa[1]))
		df<-2*(ncol(X)*(ncol(X)-1)/2+ncol(X))
		ii<-grep("Log likelihood without varA",temp)
		aa<-strsplit(strsplit(temp[ii],"=")[[1]][2]," ")[[1]]
		aa<-aa[aa!=""]
		nonVa.logLik<-as.numeric(sub(",","",aa[1]))
		nonVa.df<-ncol(X)*(ncol(X)-1)/2+ncol(X)		
		ChiSq<-2*(logLik-nonVa.logLik)
		P<-pchisq(ChiSq,df-nonVa.df,lower.tail=FALSE)
		if(hasArg(cleanup)) cleanup<-list(...)$cleanup
		else cleanup<-TRUE
		if(cleanup){
			system("rm infile",show.output.on.console=FALSE)
			system("rm intree",show.output.on.console=FALSE)
			system("rm outfile",show.output.on.console=FALSE)
		}
		if(!quiet) 	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		return(list(VarA=VarA,VarE=VarE,VarA.Regressions=VarA.Regressions,
			VarA.Correlations=VarA.Correlations,
			VarE.Regressions=VarE.Regressions,
			VarE.Correlations=VarE.Correlations,
			nonVa.VarE=nonVa.VarE,
			nonVa.VarE.Regressions=nonVa.VarE.Regressions,
			nonVa.VarE.Correlations=nonVa.VarE.Correlations,
			logLik=logLik,df=df,nonVa.logLik=nonVa.logLik,
			nonVa.df=nonVa.df))
	}
}
