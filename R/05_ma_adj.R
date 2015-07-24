### adjM(): Used by intensAdj(), performs LOESS regression on ordered MA-vals
### Input: ordered dependent variable $A$ with corresponding $M$
### Returns: adjusted $M$ values, $M_t^*$
adjM<-function(ordered.M,ordered.A)
{
	MA.finites<-is.finite(ordered.M) #only include values $>-\infty$
	finites.M<-ordered.M[MA.finites]
	finites.A<-ordered.A[MA.finites]
	
	MAloess<-loess(finites.M~finites.A
		,span=0.40,degree=2,family="symmetric",normalize=FALSE)

	finites.M<-finites.M-MAloess$fitted # make adjustments 
	ordered.M[MA.finites]<-finites.M # and return adjusted values
	return(ordered.M)
}	

### intensAdj(): Perform MA adjustment on two vectors
### Input: Two spectra vectors $F_1$ and $F_2$
### Returns: MA adjusted $F_1$ and $F_2$ values, $F_1^*$ and $F_2^*$ respectively
intensAdj<-function(F1,F2)
{
	t1<-proc.time()[3] ### get start time
	V1<-log2(F1) # Will produce $-\infty$ for $\log_2(0)$
	V2<-log2(F2)
	M<-V1-V2
	A<-(V1+V2)/2
	### $A$ is the dependent regression variable,
	### ordering required for adjM function
	ordered.indxs<-order(A)
	ordered.A<-A[ordered.indxs]
	ordered.M<-M[ordered.indxs]
	ordered.M<-adjM(ordered.M,ordered.A)
	### get indexes of original ordering
	orig.order<-order(ordered.indxs)
	M.dash<-ordered.M[orig.order]
	
	orig.finites<-is.finite(M) #update values requiring updating
	F1[orig.finites]<-2^(A[orig.finites]+M.dash[orig.finites]/2)
	F2[orig.finites]<-2^(A[orig.finites]-M.dash[orig.finites]/2)
	
	delta.t<-sprintf("%.2f",proc.time()[3]-t1) ### time elapsed
	cat("Completed MA Normalisation in",delta.t,"seconds \n")
	return(list(F1adj=F1,F2adj=F2))
}	