### m_adj(): Used by ma_adj(), performs LOESS regression on ordered MA-vals
### Input: ordered dependent variable $A$ with corresponding $M$
### Returns: adjusted $M$ values, $M_t^*$
m_adj<-function(ordered_M,ordered_A)
{
	MA_finites<-is.finite(ordered_M) #only include values $>-\infty$
	finites_M<-ordered_M[MA_finites]
	finites_A<-ordered_A[MA_finites]
	
	MAloess<-loess(finites_M~finites_A
		,span=0.40,degree=2,family="symmetric",normalize=FALSE)

	finites_M<-finites_M-MAloess$fitted # make adjustments 
	ordered_M[MA_finites]<-finites_M # and return adjusted values
	return(ordered_M)
}

### ma_adj(): Perform MA adjustment on two vectors
### Input: Two spectra vectors $F_1$ and $F_2$
### Returns: MA adjusted $F_1$ and $F_2$ values, $F_1^*$ and $F_2^*$ respectively
ma_adj<-function(F1,F2)
{
	t1<-proc.time()[3] ### get start time
	V1<-log2(F1) # Will produce $-\infty$ for $\log_2(0)$
	V2<-log2(F2)
	M<-V1-V2
	A<-(V1+V2)/2
	### $A$ is the dependent regression variable,
	### ordering required for m_adj function
	ordered_indxs<-order(A)
	ordered_A<-A[ordered_indxs]
	ordered_M<-M[ordered_indxs]
	ordered_M<-m_adj(ordered_M,ordered_A)
	### get indexes of original ordering
	orig_order<-order(ordered_indxs)
	M_dash<-ordered_M[orig_order]
	
	orig_finites<-is.finite(M) #update values requiring updating
	F1[orig_finites]<-2^(A[orig_finites]+M_dash[orig_finites]/2)
	F2[orig_finites]<-2^(A[orig_finites]-M_dash[orig_finites]/2)
	
	delta_t<-sprintf("%.2f",proc.time()[3]-t1) ### time elapsed
	cat("Completed MA Normalisation in",delta_t,"seconds \n")
	return(list(F1adj=F1,F2adj=F2))
}