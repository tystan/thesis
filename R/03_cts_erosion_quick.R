### x=$X$ are unevenly (or evenly) spaced locations of the intensities
### f=$f$ are the corresponding intensity values
### k is the SE length (size of $B$)
erode_cts_quick<-function(x,f,k){
	nx<-length(x)
	x_span<-x[nx]-x[1]
	t1<-proc.time()[3]
	isAppend<-FALSE
	if(k>=x_span){
		cat("Warning: structuring element spans the entire input set \n")
		cat("The input f vector has been output \n")
		return(f)
	}else{
		# check whether series is a length that is a multiple of k
		# if not, add a single point ($x_1+mk$,$\infty$) to series for algorithm,
		# then remove point at end
		m<-ceiling(x_span/k)
		mk<-m*k
		if(!((x[1]+mk) == x[nx])){
			x<-c(x,x[1]+mk)
			f<-c(f,+Inf)
			isAppend<-TRUE
			nx<-nx+1
			x_span<-x[nx]-x[1]
		}
		# create $\Theta$
		k_blocks<-c(0,findInterval(x,seq(x[1],x[1]+(m-1)*k,by=k)),m+1)
		# p and q are used as current $\theta$ values running along the $\Theta$ vector
		p<-k_blocks[1]
		q<-k_blocks[nx+2]
		# intialise $g,h$
		g<-rep(0,nx)
		h<-rep(0,nx)
		r_min<-rep(0,nx)
		i<-1
		j<-nx
		while(i<=nx){
			this_p<-k_blocks[i+1]
			this_q<-k_blocks[j+1]
			g[i]<-ifelse(p==this_p,min(g[i-1],f[i]),f[i])
			h[j]<-ifelse(q==this_q,min(h[j+1],f[j]),f[j])
			p<-this_p
			q<-this_q
			i<-i+1
			j<-j-1
		}
		k0<-k/2
		# fast way to determine upper and lower indexes in R to avoid looping
		lo_bounds<-nx-rev(findInterval(rev(-x),rev(-(x+k0))))+1
		hi_bounds<-findInterval(x+k0,x)
		# case 3
		r_min<-pmin(h[lo_bounds],g[hi_bounds])
		# case 1: $\theta_{w_{i}^{\triangledown}}=\theta_{w_{i}^{\vartriangle}}+1$
		which_lo<-which(k_blocks[lo_bounds]==k_blocks[hi_bounds+1])
		r_min[which_lo]<-h[lo_bounds[which_lo]]
		# case 2: $\theta_{w_{i}^{\triangledown}}+1=\theta_{w_{i}^{\vartriangle}}$
		which_hi<-which(k_blocks[lo_bounds+1]==k_blocks[hi_bounds+2])
		r_min[which_hi]<-g[hi_bounds[which_hi]]
		if(isAppend) r_min<-r_min[-nx]
		cat("Completed morphological erosion (cts scale) in"
			,sprintf("%.2f",proc.time()[3]-t1),"seconds \n")
		return(r_min)
	}
}