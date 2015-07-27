### f=$f$ are evenly spaced intensity values
### k is the SE length (size of $B$)
erode<-function(f,k) #f are the intensities
{ 
   if(!(k %% 2)) return(NULL) #SE must be of odd length
   nx<-length(f)
   erode<-rep(0,nx)
   # k0 is window coverage to the left and right of centre
   k0<-(k-1)/2
   for(i in 1:nx) #for each m/z point across the spectrum
      erode[i]<-min(f[max(1,i-k0):min(nx,i+k0)])
   return(erode)
} 
### dilate is the same as erode except use max instead of min, or:
dilate<-function(f,k) return(-erode(-f,k))
### as defined $\tau_B\left(f\right)=f-\left(f\ominus{B}\right)\oplus{B}$
tophat<-function(f,k) return(f-dilate(erode(f,k),k))