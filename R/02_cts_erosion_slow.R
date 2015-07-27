### get_lo_bounds() gets the index values of x[i]-k/2 for x[i], $i=1,2,...,n_x$.
### returns a vector of indexes length $n_x$
### used by erode_cts_slow()
get_lo_bounds<-function(x,k)
{
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
### same as get_lo_bounds() but for x[i]+k/2
get_hi_bounds<-function(x,k)
{
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
### x=$X$ are unevenly (or evenly) spaced locations of the intensities
### f=$f$ are the corresponding intensity values
### k is the SE length (size of $B$)
erode_cts_slow<-function(x,f,k)
{
	nx<-length(x)
	r_min<-rep(0,nx)
	# slow way to get LO and HI, removed to make comparison with erode_cts_quick() fair
	#LO<-get_lo_bounds(x,k)
	#HI<-get_hi_bounds(x,k)
	k0<-k/2
	# fast way to determine upper and lower indexes in R to avoid looping
	LO<-nx-rev(findInterval(rev(-x),rev(-(x+k0))))+1
	HI<-findInterval(x+k0,x)
	for(i in 1:nx) r_min[i]<-min(f[LO[i]:HI[i]])
	return(r_min)
}
