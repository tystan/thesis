erode<-function(f,sesize) #f are the intensities
{ 
   if(!(sesize %% 2)) return(NULL) #SE must be of odd length
   nx<-length(f); erode<-rep(0,nx); halfse<-(sesize-1)/2;
   for(i in 1:nx) #for each m/z point across the spectrum
      erode[i]<-min(f[max(1,i-halfse):min(nx,i+halfse)])
   return(erode)
} 
### dilate is the same as erode except use max instead of min ...OR...
dilate<-function(f,sesize) return(-erode(-f,sesize))
### as defined $\tau_B\left(f\right)=f-\left(f\ominus{B}\right)\oplus{B}$
tophat<-function(f,sesize) return(f-dilate(erode(f,sesize),sesize))