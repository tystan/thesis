
get.lo.bounds<-function(x,k){
	nx <- length(x)
	k0 <- k/2
	LO <- rep(0,nx)
	i <- 2
	i_lo <- 1
	LO[i_lo] <- i_lo # x[1] has lower bound 1
	while(i<=nx)
	{
		if((x[i]-x[i_lo])<=k0){
			LO[i] <- i_lo
			i <- i+1
		}else{
			i_lo <- i_lo+1 
		}
	}
	return(LO)
}
get.hi.bounds<-function(x,k){
	nx <- length(x)
	k0 <- k/2
	HI <- rep(0,nx)
	i <- nx-1
	i_hi <- nx
	HI[i_hi] <- i_hi # x[nx] has upper bound nx
	while(i>0)
	{
		if((x[i_hi]-x[i])<=k0){
			HI[i] <- i_hi
			i <- i-1
		}else{
			i_hi <- i_hi-1 
		}
	}
	return(HI)
}
erode.cts.slow<-function(x,f,k){
	nx<-length(x)
	r_min<-rep(0,nx)
	LO<-get.lo.bounds(x,k)
	HI<-get.hi.bounds(x,k)
	for(i in 1:nx) r_min[i]<-min(f[LO[i]:HI[i]])
	return(r_min)
}
