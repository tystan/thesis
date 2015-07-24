erode.quick<-function(f,se.size){
	nx<-length(f)
	k<-se.size
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
		k.left<-(k-1)/2
		add.pix<-k-(nx%%k)
		isMultiple<-(add.pix==k)
		if(!isMultiple){
			f<-c(f,rep(+Inf,add.pix))
			rem.indxs<-(nx+1):(nx+add.pix)
			nx<-nx+add.pix
		}
		g<-rep(0,nx)
		h<-rep(0,nx)
		r<-rep(0,nx)
		j<-nx
		for(i in 1:nx){
			if(i%%k==1){
				g[i]<-f[i]
			}else{
				g[i]<-min(g[i-1],f[i])
			}
			if(j%%k==0){
				h[j]<-f[j]
			}else{
				h[j]<-min(h[j+1],f[j])
			}
			j<-j-1
		}
		r[1:(k.left+1)]<-g[(k.left+1):k]
		r[(k.left+2):(nx-k.left-1)]<-pmin(g[(k+1):(nx-1)],h[2:(nx-k)])
		r[(nx-k.left):nx]<-h[(nx-k+1):(nx-k.left)]
		if(!isMultiple) r<-r[-rem.indxs]
		
		delta.t<-sprintf("%.2f",proc.time()[3]-t1) ### time elapsed
		cat("Completed morphological erosion in",delta.t,"seconds\n")
		return(r)
	}
}