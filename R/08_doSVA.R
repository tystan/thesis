#################################### FUNCTION ################################
#### doSVA: Perform SVA using the model:
####        $Y_j = \mu_j + X\alpha_j + Z\beta_j + W\delta_j + \mathbf{e}_j$
#################################### INPUTS ##################################
#### Y: is a $n \times p$ matrix, where each p columns are regressed
#### Intecept: boolean; do we want to fit a mean value? (yes, in most cases)
#### X: is a $n \times d_{\alpha}$ design matrix of the factors of interest
#### Z: is a $n \times d_{\beta}$ design matrix of the incidental experimental factors
#### nosigsv: the number (referred to as $H$ in some papers) of significant 
####         eigen vecs if $NULL$, the function will determine. If less than 
####         $1$, no $W$ computed
#### verbose: boolean, whether the surragate variable matrix, $W$ is returned
#### seed: an integer to feed into 'set.seed()' for reproducable results
#################################### OUTPUTS #################################
#### Ytilde: the Y matrix with $Z\beta_j + W\delta_j$ removed
#### pvals: the p-values of Ytilde regressed on $\mu_j + X\alpha_j + Z\beta_j + W\delta_j$
#### tvals: the corresponding t-statistics
#### betas: the corresponding $\alpha_j,\beta_j,\delta_j$ estimates
#### paramlabels: a combination of I (intercept), X, Z, W to signify the 
####         relevent rows of p-vals/tvals/betas
#### W: the eigen vectors matrix
#### H: the number of columns of W (used eigen-vectors)
###############################################################################
doSVA<-function(
	Y,Intercept=TRUE,X=NULL,Z=NULL,nosigsv=NULL,verbose=FALSE,seed=NULL
){
	n<-nrow(Y)
	thisInt<-IXZ<-NULL
	if(Intercept) thisInt<-matrix(1,nrow=n,ncol=1,dimnames=list(NULL,"Intcpt"))
	if(is.null(thisInt) && is.null(X) && is.null(Z))
	{
		cat("At least one of: Intercept, X and Z must be specified \n")
		return(NULL)
	} else IXZ<-cbind(thisInt,X,Z)
	kparam<-ncol(IXZ)
	colmarkers<-rep("",kparam)
	indx<-0
	if(!is.null(thisInt)) colmarkers[indx<-indx+1]<-"I"
	if(!is.null(X)) colmarkers[(indx<-indx+1):(indx<-indx+ncol(X)-1)]<-"X"
	if(!is.null(Z)) colmarkers[(indx<-indx+1):(indx<-indx+ncol(Z)-1)]<-"Z"

	RIXZ<-multReg(Y,IXZ,createNAvals=TRUE,seed=seed) 
	thissvd<-svd(RIXZ$RES)
	
	W<-H<-NULL
	if(is.null(nosigsv)){
		H<-getH(RIXZ$RES,IXZ,nullsig=0.1,verbose=FALSE)
		if(H<1) cat("No significant surrogate variables found \n")
	}else{
		H<-nosigsv
	}
	if(H<1){
		cat("No surrogate variables will be used \n")
	}else{
		cat("Using H=",H," significant surrogate variables \n",sep="")
		W<-as.matrix(thissvd$u[,1:H]) 
		colnames(W)<-paste("W",1:H,sep="")
		colmarkers<-c(colmarkers,rep("W",H))
	}
	IXZW<-cbind(IXZ,W)
	Rtilde<-multReg(Y,IXZW)
	removecols<-colmarkers %in% c("Z","W")
	ZBetaWDelta<-0
	if(sum(removecols)) ZBetaWDelta<-as.matrix(IXZW[,removecols]) %*% 
								as.matrix(Rtilde$BETA[removecols,])
	Ytilde<-Y-ZBetaWDelta
	if(verbose) return(list(Ytilde=Ytilde,paramlabels=colmarkers,W=W,H=H)) 
	else return(Ytilde)
}