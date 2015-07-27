### f=$f$ are evenly spaced intensity values
### k is the SE length (size of $B$)
erode_quick<-function(f,k){
	nx<-length(f)
	t1<-proc.time()[3] ### get start time
	if(k>=nx){
		cat("Warning: structuring element is >= in length as the input\n")
		cat("The input vector has been output \n")
		return(f)
	}else{
		if((k%%2) != 1){
			k<-k-1
			cat("Structring Element not symmetric, using SE length -1 =",k,"\n")
		}
		# k0 is window coverage to the left and right of centre
		k0<-(k-1)/2
		# check whether series is a length that is a multiple of k
		# if not, add points to series for algorithm then remove at end
		add.pix<-k-(nx%%k)
		isMultiple<-(add.pix==k)
		if(!isMultiple){
			f<-c(f,rep(+Inf,add.pix))
			rem.indxs<-(nx+1):(nx+add.pix)
			nx<-nx+add.pix
		}
		# intialise $g,h$
		g<-rep(0,nx); h<-rep(0,nx);
		r_min<-rep(0,nx)
		j<-nx
		# compute $g,h$ values
		for(i in 1:nx){
			g[i]<-ifelse(i%%k==1,f[i],min(g[i-1],f[i]))
			h[j]<-ifelse(j%%k==0,f[j],min(h[j+1],f[j]))
			j<-j-1
		}
		# only $g$ values are required at the left
		r_min[1:(k0+1)]<-g[(k0+1):k]
		# vectorised min calculations
		r_min[(k0+2):(nx-k0-1)]<-pmin(g[(k+1):(nx-1)],h[2:(nx-k)])
		# only $h$ values are required at the left
		r_min[(nx-k0):nx]<-h[(nx-k+1):(nx-k0)]
		if(!isMultiple) r_min<-r_min[-rem.indxs]
		delta.t<-sprintf("%.2f",proc.time()[3]-t1) ### time elapsed
		cat("Completed morphological erosion in",delta.t,"seconds\n")
		return(r_min)
	}
}