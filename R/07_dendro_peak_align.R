###############################################################################
##################################### Function ################################
###############################################################################
# 
# dendro_peak_align(): for peak list data, create successive N- and M-alignments
#              until all spectra are aligned.
# 
###############################################################################
####################################### Input #################################
###############################################################################
#
# msD: MS Data, a $T x n$ matrix of MS intensities. One column per spectra.
# peaklistlist: see below
# in.param: [ D expon lambda G maxM ]-tuple as a vector
#
###############################################################################
###################################### Output #################################
###############################################################################
#
### A list containing the following elements:
# dendro: 
# stepwise.peaks: a list where each element is the successive amalgamation data
# amalpeaks: the final matrix of aligned peaks. Columns are named spectra, rows
#    are
#
###############################################################################

### peaklistlist=
# [[1]]
# [t_{1,1} t_{1,2} ... t_{1,n_1} ]
# [x_{1,1} x_{1,2} ... x_{1,n_1} ]
# 
# [[2]]
# [t_{2,1} t_{2,2} ... t_{2,n_2} ]
# [x_{2,1} x_{2,2} ... x_{2,n_2} ]
# 
# .
# .
# .
# 
# [[N]]
# [t_{N,1} t_{N,2} ... t_{N,n_N} ]
# [x_{N,1} x_{N,2} ... x_{N,n_N} ]
# 
# where t_{i,j} is the time point j-th peak for the i-th spectrum
# where x_{i,j} is the vector of intensities (nComp long i.e. nComp x 1 matrix) 
# 				for the j-th peak for the i-th spectrum
# NB: each list item is a ($n_{Comp}$+1) x $n_N$ matrix

dendro_peak_align<-function(msD,peaklistlist,in.param)
{

	D<-in.param[1]
	nC<-nrow(peaklistlist[[1]])-1
	expon<-in.param[2]
	lambda<-in.param[3]
	G<-in.param[4]
	maxM<-in.param[5]
	
	nPat<-length(peaklistlist)
	Pats<-1:nPat
	
	cat("Calculating merge sequence for spectra \n")
		
	fordist<-t(msD$intensity)
	hc<-hclust(as.dist(fordist,diag=FALSE,upper=FALSE),"average")
	
	### find amalgamation sequence
	### see ?hclust for information on the merge matrix:
	#    "an n-1 by 2 matrix. Row i of merge describes the merging of clusters 
	#    at step i of the clustering. If an element j in the row is negative, 
	#    then observation -j was merged at this stage. If j is positive then 
	#    the merge was with the cluster formed at the (earlier) stage j of the 
	#    algorithm. Thus negative entries in merge indicate agglomerations of 
	#    singletons, and positive entries indicate agglomerations of 
	#    non-singletons."
	amalg<-hc$merge

	nAmal<-nrow(amalg)
	# alignment of peaklistlist (a.pll)
	a.pll<-vector(length=nAmal,mode="list")
	Npeaks<-NULL
	Npeaklist<-NULL
	Mpeaks<-NULL
	Mpeaklist<-NULL
	
	# start
	for(aindx in 1:nAmal)
	{
	
		# if patsToGetN or patsToGetM are positive,
		#    ... it is a single spectrum (1-alignment)
		# if negative, it is a previous N/M-alignment (N,M>1) 
		patsToGetN<--amalg[aindx,1] 
		patsToGetM<--amalg[aindx,2]
		
		printPatsN<-sprintf("%03d",patsToGetN)
		printPatsM<-sprintf("%03d",patsToGetM)
		amalg.str<-"Amalgamting patient"
		if(patsToGetN>0 && patsToGetM>0){
			cat(amalg.str,"s ",printPatsN," and ",printPatsM,"\n",sep="")
		}else if(patsToGetN>0 && patsToGetM<0){ 
			cat(amalg.str,printPatsN,"to previously amalgamated patients\n")
		}else if(patsToGetN<0 && patsToGetM>0){ 
			cat(amalg.str,printPatsM,"to previously amalgamated patients\n")
		}else{ 
			cat("Amalgamting two clusters of previously amalgamated patients\n")
		}
		
		### prepare N-Alignment data
		if(patsToGetN>0){ # if a single spectrum (1-alignment)
			Npeaks<-matrix(1:ncol(peaklistlist[[patsToGetN]]),ncol=1)
			Npeaklist<-peaklistlist[patsToGetN]
		}else{ # if a previously aligned N-alignment (N>1)
			Npeaks<-a.pll[[-patsToGetN]]
			patsToGetN<-as.numeric(colnames(Npeaks))
			Npeaklist<-peaklistlist[patsToGetN]
		}

		### prepare M-Alignment data
		if(patsToGetM>0){ # if a single spectrum (1-alignment)
			Mpeaks<-matrix(1:ncol(peaklistlist[[patsToGetM]]),ncol=1)
			Mpeaklist<-peaklistlist[patsToGetM]
		}else{ # if a previously aligned M-alignment (M>1)
			Mpeaks<-a.pll[[-patsToGetM]]
			patsToGetM<-as.numeric(colnames(Mpeaks))
			Mpeaklist<-peaklistlist[patsToGetM]
		}
		
		### use Wmatrix() function
		Wm<-Wmatrix(Npeaks,Npeaklist,Mpeaks,Mpeaklist,D,expon,lambda)
		### use S-W alignment function to estimate maximum path
		### see: https://code.google.com/p/swalign/
		estPM<-SWalign(Wm,G,maxM) 
		### estPM is a data.frame of (i,j) locations of the maximum path
		### the data.frame is 2 columns for i,j points
		
		nN<-ncol(Npeaks) ### no. of peaks in N-align
		nM<-ncol(Mpeaks) ### no. of peaks in M-align
		nK<-nrow(estPM)  ### no. of peaks in new N:M-align
		### apllTemp: 
		###      (a)lignment of (p)eak (l)ist (l)ist, (temp)orary
		### Matrix of peak indicators. The $n_K$ rows represent the $n_K$ peaks 
		###    from the N:M-alignment. 
		### Entries apllTemp[i,j] are ==
		### { 0 if that N:M-aligned peak does not exist in spec $j$ (column $j$)
		### { _else_ a non-zero indicator, the peak number from within the  
		###                           1-alignment from spectrum $j$ (column $j$)
		apllTemp<-matrix(0,nrow=nK,ncol=nN+nM)
		mzValsTemp<-NULL
		AveMzValsTemp<-NULL
		for(n.k in 1:nK)
		{
			if(estPM[n.k,1]>0) # if the peak exists in the N-alignment
			{
				# transfer peak info from N-align to new N:M-align matrix
				apllTemp[n.k,1:nN]<-Npeaks[estPM[n.k,1],]
				for(i in 1:nN) if(apllTemp[n.k,i]>0) mzValsTemp<-
						c(mzValsTemp,Npeaklist[[i]][1,apllTemp[n.k,i]])
			}
			if(estPM[n.k,2]>0) # if the peak exists in the M-alignment
			{
				# transfer peak info from N-align to new N:M-align matrix
				apllTemp[n.k,(nN+1):(nN+nM)]<-Mpeaks[estPM[n.k,2],]
				for(i in (nN+1):(nN+nM)) if(apllTemp[n.k,i]>0) mzValsTemp<-
						c(mzValsTemp,Mpeaklist[[i-nN]][1,apllTemp[n.k,i]])
			}
			# get ave m/z of all aligned peaks
			AveMzValsTemp<-c(AveMzValsTemp,mean(mzValsTemp)) 
			mzValsTemp<-NULL
		}
		### change row order if averaging m/z has changed peak location order
		mzReOrder<-order(AveMzValsTemp)
		apllTemp<-apllTemp[mzReOrder,]
		
		allPat<-c(patsToGetN,patsToGetM)
		colnames(apllTemp)<-allPat
		### clean up N:M-alignment to preserve spectrum order
		patOrder<-order(allPat)
		apllTemp<-apllTemp[,patOrder]
		
		a.pll[[aindx]]<-apllTemp
	}
	### return list() object of peak amalgamation/alignment,
	###     including intermediate steps
	outlist<-list(dendro=hc,stepwise.peaks=a.pll[-nAmal],amalpeaks=a.pll[[nAmal]])
	return(outlist)

}