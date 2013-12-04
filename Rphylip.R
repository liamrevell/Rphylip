## function to crop to first n characters a vector of strings
## written by Liam J. Revell 2013

crop<-function(x,n=1) sapply(x,function(x) strsplit(x,"")[[1]][1:n])

## function to change vector to integers
## written by Liam J. Revell 2013

to.integers<-function(x){
	types<-sort(unique(x))
	ii<-as.integer(1:length(types)-1)
	ii[sapply(x,function(x,y) which(y==x),y=types)]
}

## calls threshml from PHYLIP (Felsenstein 2013)
## written by Liam J. Revell 2013

Rthreshml<-function(tree,X,types=NULL,path=NULL,...){
	if(is.null(path)) path<-findPath("threshml")
	if(is.null(path)) stop("No path provided and was not able to find path to threshml")
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	if(is.null(types)){
		types<-sapply(X,class)
		types[types=="numeric"]<-"continuous"
		types[types%in%c("factor","character")]<-"discrete"
	}
	types<-crop(types)
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile"))==0) return(NULL)
	tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y)[1],y=rownames(X))
	## this is for the current idiosyncratic tree input file requirement of threshml
	text<-write.tree(tree)
	text<-strsplit(text,"")[[1]]
	text<-paste(paste(text[1:(length(text)-1)],collapse=""),"0.00000000;\n",sep="")
	write(text,file="intree")
	if(any(types=="c")) write.continuous(X[,crop(types)=="c"])
	if(any(types=="d")) write.dna(apply(X[,crop(types)=="d"],2,to.integers),append=any(types=="c"))
	## start populating arguments
	oo<-c("r")
	if(!any(types=="d")) oo<-c(oo,"d")
	if(any(types=="c")) oo<-c(oo,"c")
	if(hasArg(burnin)){
		burnin<-list(...)$burnin
		oo<-c(oo,"b",burnin)
	}
	if(hasArg(nchains)){
		nchains<-list(...)$nchains
		oo<-c(oo,"n",nchains)
	}
	if(hasArg(ngen)){
		ngen<-list(...)$ngen
		oo<-c(oo,"s",ngen)
	}
	if(hasArg(proposal)){
		proposal<-list(...)$proposal
		oo<-c(oo,"p",proposal)
	}
	if(hasArg(lrtest)) lrtest<-list(...)$lrtest
	else lrtest<-TRUE
	if(lrtest) oo<-c(oo,"t")
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y",sample(seq(1,99999,by=2),1))
	system("touch outfile")
	system(paste(path,"/threshml",sep=""),input=oo)
	temp<-readLines("outfile")
	
##		ii<-grep("Contrasts",temp)
##		Contrasts<-matrix(NA,tree$Nnode,ncol(X))
##		for(i in 1:tree$Nnode){
##			x<-strsplit(temp[i+ii+2]," ")[[1]]
##			Contrasts[i,]<-as.numeric(x[x!=""])
##		}
##		ii<-grep("Covariance",temp)
##		Covariance_matrix<-matrix(NA,ncol(X),ncol(X))
##		for(i in 1:ncol(X)){
##			x<-strsplit(temp[i+ii+2]," ")[[1]]
##			Covariance_matrix[i,]<-as.numeric(x[x!=""])
##		}
##		ii<-grep("Regressions",temp)
##		Regressions<-matrix(NA,ncol(X),ncol(X))
##		for(i in 1:ncol(X)){
##			x<-strsplit(temp[i+ii+2]," ")[[1]]
##			Regressions[i,]<-as.numeric(x[x!=""])
##		}
##		ii<-grep("Correlations",temp)
##		Correlations<-matrix(NA,ncol(X),ncol(X))
##		for(i in 1:ncol(X)){
##			x<-strsplit(temp[i+ii+2]," ")[[1]]
##			Correlations[i,]<-as.numeric(x[x!=""])
##		}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup) cleanFiles(c("infile","intree","outfile"))
	if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
##		return(list(Contrasts=Contrasts,Covariance_matrix=Covariance_matrix,
##			Regressions=Regressions,Correlations=Correlations))
	return(NULL)
}

## function to write continuous characters to file
## written by Liam J. Revell 2013

write.continuous<-function(X,append=FALSE){
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile",append=append)
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## calls dnamlk from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnamlk<-function(X,path=NULL,...){
	Rdnaml(X,path,clock=TRUE,...)
}	

## clean up files
## written by Liam J. Revell 2013

cleanFiles<-function(fs){
	if(.Platform$OS.type=="windows") for(i in 1:length(fs)) system(paste("rm",fs[i],sep=" "),show.output.on.console=FALSE)
	else for(i in 1:length(fs)) system(paste("rm",fs[i],sep=" "))
}

## sets up PHYLIP in Mac OS X (based on http://evolution.gs.washington.edu/phylip/install.html)
## written by Liam J. Revell

setupOSX<-function(path=NULL){
	if(.Platform$OS.type!="unix") stop("this function is for Mac OS X only")
	if(is.null(path)){
		## check /Applications for path to PHYLIP
		ll<-list.files("/Applications/")
		ii<-grep("phylip",ll)
		if(length(ii)>0) path<-paste("/Applications/",ll[ii],sep="")
		else stop("was not able to find path to phylip installation")
	}
	if(strsplit(path,"")[length(strsplit(path,""))]=="/"){
		path<-strsplit(path,"")
		path<-paste(path[2:length(path)-1],collapse="")
	} 
	system(paste("cp ",path,"/src/linkmac ",path,"/exe/linkmac",sep=""))
	system(paste("chmod +x ",path,"/exe/linkmac",sep=""))
	system(paste("cd ",path,"/exe/\n./linkmac",sep=""))
}		

## calls neighbor from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rneighbor<-function(D,path=NULL,...){
	if(class(D)=="dist"||class(D)=="data.frame") D<-as.matrix(D)
	D<-D[rownames(D),rownames(D)]
	if(is.null(path)) path<-findPath("neighbor")
	if(is.null(path)) stop("No path provided and was not able to find path to neighbor")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"nj"
	if(method=="NJ"||method=="nj") method<-"nj"
	else if(method=="UPGMA"||method=="upgma"){
		method<-"upgma"
		oo<-c(oo,"n")
	} else {
		cat("\nWarning:\n  method not recognized - using method=\"NJ\"\n")
		method="nj"
	}
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order) oo<-c(oo,"j",sample(seq(1,99999,by=2),1))
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.distances(D)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/neighbor",sep=""),input=oo)
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(D),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(D))
		cat("\n")
	}
	tree$tip.label<-rownames(D)[as.numeric(tree$tip.label)]	
	if(hasArg(outgroup)){
		if(method=="nj"){
			outgroup<-list(...)$outgroup
			tree<-outgroup.root(tree,outgroup,quiet)
		} else { 
			cat("\nWarning:\n  outgroup rooting not permitted for method = \"upgma\"\n")
			cat("  tree already rooted\n\n")
		}
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup) cleanFiles(c("infile","outfile","outtree"))
	return(tree)
}

## write distance matrix to file in PHYLIP format
## written by Liam J. Revell 2013

write.distances<-function(D){
	write(paste("    ",nrow(D),sep=""),file="infile")
	for(i in 1:nrow(D)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(D[i,],collapse=" "),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## attempt to find path to PHYLIP executable
## written by Liam J. Revell 2013

findPath<-function(string){
	if(.Platform$OS.type=="windows"){
		## first, check current directory
		ll<-list.files()
		ii<-grep(string,ll)
		if(length(ii)>0) 
			if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
				return(".")
		## check C:/Program Files
		ll<-list.files("C:/Program Files/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Program Files/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))
		}
		## check C:/Progam Files (x86)
		ll<-list.files("C:/Program Files (x86)/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Program Files (x86)/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep="")))
					return(shortPathName(dd))
		}
		## check C:/Users/Username
		uu<-strsplit(getwd(),"")[[1]]
		ii<-grep("/",uu)
		uu<-paste(uu[(ii[2]+1):(ii[3]-1)],collapse="")
		ll<-list.files(paste("C:/Users/",uu,sep=""))
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Users/",uu,"/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))
		}
		## check C:/Users/Username/Documents
		ll<-list.files(paste("C:/Users/",uu,"/Documents",sep=""))
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Users/",uu,"/Documents/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))

		}
		return(NULL)
	} else if(.Platform$OS.type=="unix"){
		## first, check current directory
		ll<-list.files()
		ii<-grep(string,ll)
		if(length(ii)>0) 
			if(any(ll[ii]==string)) 
				return(".")
		## check /Applications
		ll<-list.files("/Applications/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("/Applications/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(dd)
		}
		return(NULL)
	} else return(NULL)
}

## function writes DNAbin to file in PHYLIP format with numbers as labels
## written by Liam J. Revell 2013

write.dna<-function(X,append=FALSE){
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile",append=append)
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[i,],collapse=""),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## function to outgroup root
## written by Liam J. Revell 2013

outgroup.root<-function(tree,outgroup,quiet){
	if(class(tree)=="phylo") tree<-root(tree,outgroup)
	else if(class(tree)=="multiPhylo"){
		tree<-lapply(tree,root,outgroup=outgroup)
		class(tree)<-"multiPhylo"
	}
	if(!quiet){
		cat("Rooted tree(s) with the outgroup\n")
		cat("------------------------\n")
		cat(paste(paste(outgroup,collapse=", "),"\n\n"))
	}
	return(tree)
}

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

Rdnapars<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnapars")
	if(is.null(path)) stop("No path provided and was not able to find path to dnapars")
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
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnapars",sep=""),input=oo)
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
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
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## calls contml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontml<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("contml")
	if(is.null(path)) stop("No path provided and was not able to find path to contml")
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
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)	
	}
	tree$logLik<-logLik
	return(tree)
}

## calls dnaml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnaml<-function(X,path=NULL,...){
	if(hasArg(clock)) clock<-list(...)$clock
	else clock<-FALSE
	exe<-if(clock) "dnamlk" else "dnaml"
	if(is.null(path)) path<-findPath(exe)
	if(is.null(path)) stop(paste("No path provided and was not able to find path to",exe))
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
	if((!speedier)&&(!clock)) oo<-c(oo,"s")
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
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	temp<-system(paste(path,"/",exe,sep=""),input=oo,show.output.on.console=(!quiet))
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
		if(!clock) tree<-outgroup.root(tree,outgroup,quiet)
		else cat("\nMolecular clock trees are already rooted!\n\nIgnoring argument outgroup.\n\n")
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(rates)) files<-c(files,"rates")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	tree$logLik<-logLik
	return(tree)
}

## function to optimize parameters of Rdnaml
## written by Liam J. Revell 2013

opt.Rdnaml<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnaml")
	if(is.null(path)) stop("No path provided and was not able to find path to dnaml")
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
	## bounds for optimization
	if(hasArg(bounds)) bounds<-list(...)$bounds
	else bounds<-list()
	bb=list(kappa=c(0.01,20),gamma=c(0.01,20),bf=cbind(rep(0.01,4),rep(1,4)))
	bb[(namc<-names(bounds))]<-bounds
	par<-c(1,10,rep(0.25,4))
	fit<-optim(par,lik,X=X,tree=tree,path=path,method="L-BFGS-B",
		lower=c(bb$kappa[1],bb$gamma[1],bb$bf[,1]),
		upper=c(bb$kappa[2],bb$gamma[2],bb$bf[,2]),control=list(fnscale=-1))
	return(list(kappa=fit$par[1],gamma=fit$par[2],
		bf=fit$par[3:6]/sum(fit$par[3:6]),logLik=fit$value))
}

## function calls contrast from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontrast<-function(tree,X,path=NULL,...){
	if(is.null(path)) path<-findPath("contrast")
	if(is.null(path)) stop("No path provided and was not able to find path to contrast")
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
		write.continuous(X)
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
		if(cleanup) cleanFiles(c("infile","intree","outfile"))
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
		if(cleanup) cleanFiles(c("infile","intree","outfile"))
		if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
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
