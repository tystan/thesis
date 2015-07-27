# require animation package for kfcv() function
library(animation)
# load compiled C code (shared object)
dyn.load("dom_feat.so")

##########################################
##########################################
#
# Pairwise case of Pareto Fronts
#     obj is the $2 \times n$ matrix of the two vectors of length $n$ for  
#			the features/observations of the 2 criteria/objective functions
#     istomin is a boolean vector of whether the criteria obj$_1$, obj$_2$  
#			are to be minimised (=TRUE), respectively
#
##########################################
##########################################

#
# This function returns a vector of 'dominated' observations (Boolean, 
# length $n$ vector) FALSE=Pareto front, TRUE=dominated observation
#

dom_features_pw<-function(obj,istomin)
{
	#if to be minimised then just make negative and maximise
	obj[,istomin]<- -obj[,istomin] 
	n<-as.integer(nrow(obj))
	obj1<-as.double(obj[,1])
	obj2<-as.double(obj[,2])
	domvec<-as.integer(rep(0,n)) #output vector
	
	return(as.logical(.C("dom_feat",n,obj1,obj2,domvec)[[4]]))
}
	
##########################################
##########################################
#
# General case
#     objmatrix is $n \times m$ matrix. $n$ features/observations and 
#     $m$ criteria/objective functions istominvec is a boolean vec
#     of length $m$ to say whether the criteria are to be minimised
#
##########################################
##########################################

#
# This function returns a vector of 'dominated' observations (Boolean,
# length $n$ vector) FALSE=Pareto front, TRUE=dominated observation
# same as pairwise but the input can take more than two objective functions
#

dom_features<-function(objmatrix,istominvec)
{
	n<-as.integer(nrow(objmatrix))
	m<-ncol(objmatrix)
	objmatrix[,istominvec]<- -objmatrix[,istominvec]
	vecdomvec<-rep(1,n)
	indxs<-combn(m,2)
	nm<-ncol(indxs)
	i<-0
	while(i<nm)
	{
		i<-i+1
		# call pairwise function, take the intersection of previous 
		# dominated observations remembering the intersection(s)
		# of dominated in the same as unions(s) of Pareto fronts
		vecdomvec<-vecdomvec * .C("dom_feat"
							,n
							,as.double(objmatrix[,indxs[1,i]])
							,as.double(objmatrix[,indxs[2,i]])
							,as.integer(rep(0,n)))[[4]]
	}
	return(as.logical(vecdomvec))
}

##########################################
##########################################
#
# Sucessive Pareto Fronts
#     noFronts is the # of pareto fronts required
#
#     fn returns a vector of length $n$
#           each element is labelled the pareto front #,  
#           0 is dominated even after noFronts found
#
##########################################
##########################################

pareto_fronts<-function(noFronts,objmatrix,istominvec)
{
	objmatrix[,istominvec]<- -objmatrix[,istominvec]
	n<-as.integer(nrow(objmatrix))
	m<-ncol(objmatrix)
	pfvec<-rep(0,n) #output vector
	#once a front is found we need to set the correponding values to $\infty$ or 
	# $-\infty$ so they won't be chosen again
	#try: as.numeric(c(TRUE,FALSE,TRUE))*2-1 to see what the next line is doing
	#if Min then set 1, ifMax then set -1 (the sign of the Inf if we find front 
	#	 and have to put to a value)
	ourInfs<-min(objmatrix)-1
	allFrontsFound<-FALSE
	
	i<-0
	while(i<noFronts && !allFrontsFound) #go thru all fronts required
	{
		i<-i+1
		df<-dom_features(objmatrix,rep(FALSE,m)) #general m obj vectors function
		# pf.i are the indexs of the output vector that need to be updated 
		# with the pareto front number
		pf.i<-(!df) & (pfvec<1) 
		pfvec[pf.i]<-i
		# re-assign values were pareto front found
		objmatrix[pf.i,]<-ourInfs
		if(all(pfvec>0)) allFrontsFound<-TRUE
	}
	return(pfvec)
}

##########################################
##########################################
#
# Leave-one-out/k-fold feature ranking
#
# returns a vector of length $n$ with values $\in (0,1]$ for feature importance		
#
##########################################
##########################################

#
# Same inputs of previous functions, with folds (aka $k$-fold cross  
# validation) and reps is the number of times we re-do the cross 
# validation fold=1 or the length of the input (i.e. n) creates 
# leave-one-out cross validation

pareto_ranking<-function(objmatrix,istominvec,noFronts=20,folds=1,reps=5)
{
	objmatrix[,istominvec]<- -objmatrix[,istominvec]
	m<-ncol(objmatrix)
	n<-nrow(objmatrix)
	pfmetric<-rep(0,n) #output vector
	nfolds<-n
	if(folds>1) nfolds<-folds
	if(nfolds==n) reps<-1
	
	blocks<-kfcv(nfolds,n)
	block.nos<-rep(1:nfolds,blocks)
	
	for(r in 1:reps)
	{
		indxs<-sample(1:n) #fresh randomisation each repetition
		k.f.mat<-cbind(indxs,block.nos) # create the fold 'blocks' of data
		for(i in 1:nfolds)
		{
			rows<-k.f.mat[k.f.mat[,2]==i,1] # find the ith fold to leave out
			# call general function with ith fold removed
			calcfronts<-pareto_fronts(noFronts,objmatrix[-rows,],rep(FALSE,m))
			# which are non-dominated
			whichnondom<-calcfronts>0
			# if you are ont the first front you get 1, second=1/2, third=1/3,
			# ..., jth=1/j, else 0
			pfmetric[-rows][whichnondom]<-
					pfmetric[-rows][whichnondom]+1/calcfronts[whichnondom]
		}
	}
	#now divide by maximum posible value i.e. (nfolds-1)*reps so output in [0,1]
	pfmetric<-pfmetric/((nfolds-1)*reps) 
	return(pfmetric)
}

###########################################################################
# Below are three metrics that can possibly describe the value of      
#   variables/fetures to discriminate betwwen classes                                                             ####
###########################################################################

# minIntraClassVar(): find the minimum WITHIN class variance of the K groups
# interClassVar(): find the variance of means/centroids of the K groups
# maxInterClassDist(): possibly correlated with interClassVar, find the dist 
#							MAXIMUM between the K group's centroids/means
                         
# The rationale of the last metric is that a variable/feature that only 
# seperates two of the K classes  is undervalued by the Fisher score because 
# it may not separate the K-2 classes remaining well.	
# ... And a separation of two classes (in conjunction with other variables) 
#      is important information for the discriminant model

###########################################################################
#### ds: a data.frame or matrix (numeric values only/factors will be dealt 
####	with as integers)  class vec: a vector correspong to the class of 
####    the rows of ds           
###########################################################################

minIntraClassVar<-function(ds,class.vec){

	dsfs<-ds
	if(!is.matrix(dsfs)) dsfs<-data.matrix(dsfs)  
	
	p<-ncol(dsfs)
	K<-length(levels(class.vec))
	n.all<-length(class.vec)
	n.i<-0
	
	intraClassVar<-Inf
	mean.j<-colMeans(dsfs)
	
	for(i in 1:K){
		true.vec<-(as.integer(class.vec)==i)
		n.i<-length(which(true.vec))
		mean.class<-colMeans(as.matrix(dsfs[true.vec,],ncol=p))
		var.class<-colSums(as.matrix((dsfs[true.vec,]-rep(mean.class,each=n.i))^2,ncol=p))
		intraClassVar<-pmin(intraClassVar,var.class/(n.i-1))
	}
	return(intraClassVar)
}

interClassVar<-function(ds,class.vec){

	dsfs<-ds
	if(!is.matrix(dsfs)) dsfs<-data.matrix(dsfs)  
	
	p<-ncol(dsfs)
	K<-length(levels(class.vec))
	
	interClassVar<-0
	mean.j<-colMeans(dsfs)
	
	for(i in 1:K){
		true.vec<-(as.integer(class.vec)==i)
		n.i<-length(which(true.vec))
		mean.class<-colMeans(as.matrix(dsfs[true.vec,],ncol=p))
		var.class<-((mean.j-mean.class)^2)
		interClassVar<-interClassVar+var.class
	}
	return(interClassVar/(K-1))
}

maxInterClassDist<-function(ds,class.vec){

	dsfs<-ds
	if(!is.matrix(dsfs)) dsfs<-data.matrix(dsfs)  
	
	p<-ncol(dsfs)
	K<-length(levels(class.vec))
	
	interClassDist<--Inf
	mean.j<-colMeans(dsfs)
	
	for(i in 1:K){
		true.vec<-(as.integer(class.vec)==i)
		mean.class<-colMeans(as.matrix(dsfs[true.vec,],ncol=p))
		interClassDist<-pmax(interClassDist,abs(mean.j-mean.class))
	}
	return(interClassDist)
}