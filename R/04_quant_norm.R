### Input: msData - the spectra intensities in matrix
###        where columns are spectra $1,2,\hdots,n$
msEmpiricalQuantNorm<-function(msData)
{
	msD<-msData
	nSpec<-ncol(msD)
	nDim<-nrow(msD)

	orders<-apply(msD, 2, order)
	reorders<-apply(orders, 2, order)
	# order each column into ascending order
	for(i in 1:nSpec) msD[,i]<-msD[orders[,i],i]
	#replace ordered columns with row means
	rmeans<-rowMeans(msD)
	for(i in 1:nSpec) msD[,i]<-rmeans
	#put back into the original order (with changed values)
	for(i in 1:nSpec) msD[,i]<-msD[reorders[,i],i]

	return(msD)
}