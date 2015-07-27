###############################################################################
##################################### Function ################################
###############################################################################
# 
# w_matrix(): Create a peak similarity matrix between an N- and M-alignment
# 
###############################################################################
####################################### Input #################################
###############################################################################
#
# Nmatchedpeaks: K x N matrix of matched pairs of the N alignment
# Npeaklists: a list object containing N matricies of  
#         (time,intensityVector) pairs: ($n_a$ x ($n_{Comp}$+1) matrix, a=1,...,N)
# Mmatchedpeaks: L x M matrix of matched pairs of the M alignment
# Mpeaklists: a list object containing M matricies of  
#         (time,intensityVector) pairs: ($n_b$ x ($n_{Comp}$+1) matrix, b=1,...,M)
#
###############################################################################
# 
# e.g. Nmatchedpeaks =
# [ 1 0 1 0 ]
# [ 0 1 2 1 ]
# [ 0 0 0 2 ]
# [ 2 2 0 0 ]
# [ 3 3 3 3 ]
# [ 0 4 0 4 ]
# [ . . . . ]
# [ . . . . ]
# [ . . . . ]
# 
# here K x N (N=4) matrix
# 
# Npeaklist=
# [[1]]
# [$t_{1,1}$ $t_{1,2}$ ... $t_{1,n_1}$ ]
# [$x_{1,1}$ $x_{1,2}$ ... $x_{1,n_1}$ ]
# 
# [[2]]
# [$t_{2,1}$ $t_{2,2}$ ... $t_{2,n_2}$ ]
# [$x_{2,1}$ $x_{2,2}$ ... $x_{2,n_2}$ ]
# 
# ...
# 
# [[N]]
# [$t_{N,1}$ $t_{N,2}$ ... $t_{N,n_N}$ ]
# [$x_{N,1}$ $x_{N,2}$ ... $x_{N,n_N}$ ]
# 
# where $t_{i,j}$ is the time point j-th peak for the i-th spectrum
# where $x_{i,j}$ is the vector of intensities ($n_{Comp}$ long)  
#      i.e. $n_{Comp}$ x 1 matrix for the j-th peak for the i-th spectrum
# NB: each list item is a $(n_{Comp}+1)$ x $n_N$ matrix
#
# the _i_ co-ord is the row in the Nmatchedpeaks
# the _a_ co-ord is the column number of Nmatchedpeaks (the spectrum number)
# the _p_ co-ord is the peak number for the _a_th spectrum
#
# the _j_,_b_ and _q_ co-ords are defined similarly for the M-alignment

w_matrix<-function(Nmatchedpeaks,Npeaklists,Mmatchedpeaks,Mpeaklists
					,D,expon,lambda){

   K<-nrow(Nmatchedpeaks)
   N<-ncol(Nmatchedpeaks)
   L<-nrow(Mmatchedpeaks)
   M<-ncol(Mmatchedpeaks)
   W<-matrix(0,nrow=K,ncol=L)
   
   for(i in 1:K){
      for(j in 1:L){
         numerator<-0
         denominator<-0
         for(a in 1:N){
            p<-Nmatchedpeaks[i,a]
            if(p>0){
               for(b in 1:M){
                  q<-Mmatchedpeaks[j,b]
                  if(q>0){
                     t_a<-Npeaklists[[a]][1,p]
                     p_a<-Npeaklists[[a]][-1,p]
                     t_b<-Mpeaklists[[b]][1,q]
                     p_b<-Mpeaklists[[b]][-1,q]
                     numerator<-numerator+
                     	PeakSim(p_a,t_a,p_b,t_b,D,expon,lambda)
                     denominator<-denominator+1
                  }
               }
            }
         }
         if(denominator>0) W[i,j]<-numerator/denominator
         else W[i,j]<-0
      }
   }
   return(W)
}