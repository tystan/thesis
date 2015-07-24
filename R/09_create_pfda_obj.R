####### FUNCTION: createPFldaobj()
### estimate parameters of PF-DA model, so that a discrim function created

####### input:
### X: a n x p matrix, of n obs and p variables
### Xclass: a vector of length n of the classes (must be a factor variable)
### priors: a vector of length K (#classes) with elements in (0,1)

createPFldaobj<-function(X,Xclass,lambdar=1,priors=NULL,alph=NULL,wts=NULL)
{
  N<-length(Xclass)
  P<-ncol(X)
  nks<-table(Xclass)
  classnames<-levels(Xclass)
  K<-length(classnames)
  
  ### if not supplied, make $\hat{\pi}_k$ data proportions
  if(is.null(priors)) priors<-nks/N
  
  if(length(priors)!=K){
    cat("The length of priors and the total number 
        of groups must be equal \n")
    return(NULL)
  }else if(is.null(alph) & (N<P)){
    cat("Alpha is suggested for n<p data \n")
  }else if(N!=nrow(X)){
    cat("The length of Xclass and the number 
        of rows in X must agree \n")
    return(NULL)
  }else if(!all(nks>1)){
    cat("There needs to be at least two obs in each 
        group for variances to be computed \n")
    return(NULL)
  }

  Xclassint<-as.integer(Xclass)
  transMeans<-colMeans(X)
  X<-X-matrix(rep(transMeans,N),nrow=N,byrow=TRUE)

  ##### create $\mu_k=\left[\mu_{k1},\hdots,\mu_{kp}\right]$ vectors
  ##### place on top of each other to get KxP matrix
  MuMat<-matrix(0,nrow=K,ncol=P)
  for(k in 1:K) MuMat[k,]<-colMeans(X[Xclassint==k,])
  MuIter<-MuMat

  ##### create $\Sigma$
  Sigma<-matrix(0,nrow=P,ncol=P)
  for(k in 1:K) 
  {
    rowuse<-which(Xclassint==k)
    Sigma<-Sigma+length(rowuse)*cov.wt(X[rowuse,],cor=FALSE
                                      ,center=TRUE,method="ML")$cov
  }
  Sigma<-Sigma/N
  ##### extract Diag elements
  sigmasqs<-diag(Sigma)
  
  if(!(is.null(alph) | is.null(wts))) sigmasqs<-sigmasqs+alph*wts
  else if(!is.null(alph)) sigmasqs<-sigmasqs+rep(alph,P)
  
  #### Now start iterative estimation of the $\ell_1$ penalised means
  #### Note "squig" is used for the ML estimates
  G<-matrix(0,nrow=K,ncol=K)
  deltatol<-1e-10
  deltaMu<-1
  maxIter<-500
  itcount<-0
  while(deltaMu>(1e-5) && itcount<maxIter)
  {
    deltaMuNumer<-0
    deltaMuDenom<-0
    itcount<-itcount+1
    
    ### For each of the features
    for(j in 1:P) 
    {
      beta.t.j<-MuIter[,j]
      musqig.j<-MuMat[,j]
      sqigY<-X[,j]
      sqigX<-matrix(0,nrow=N,ncol=K)
      G<-matrix(0,nrow=K,ncol=K)
      
      ### $\sum_{k=1}^{K-1}\sum_{k_{dash}=k+1}^K$
      for(k in 1:(K-1))
      {
        for(kdash in (k+1):K)
        {
          PFweight<-1/abs(musqig.j[k]-musqig.j[kdash])
          ### assign updated iterations, or tol value if "zero"
          muDiffIter<-max(abs(beta.t.j[k]-beta.t.j[kdash]),deltatol)
          G[k,kdash]<-G[kdash,k]<- -PFweight/muDiffIter
        }
      }
      for(k in 1:K) 
      {
        sqigX[which(Xclassint==k),k]<-1  
        ### note the diag elements of G can be calculated as the sum of the column
        G[k,k]<- -sum(G[,k])
      }
      #### $\hat{M}=\left(B^TB+\lambda\sigma_j^2G\right)^{-1}B^TJ$
      MuIter[,j]<-solve(t(sqigX)%*%sqigX+lambdar*sigmasqs[j]*G)%*%(t(sqigX)%*%sqigY)
      deltaMuNumer<-deltaMuNumer+sum(abs(MuIter[,j]-beta.t.j))
      deltaMuDenom<-deltaMuDenom+sum(abs(beta.t.j))
    }
    ### our break loop value
    deltaMu<-deltaMuNumer/deltaMuDenom
  
  }
  cat("Iterations performed to aquire a solution:",itcount
       ,"| Final tol val:",deltaMu," \n")
  
  MuIter[which(MuIter<deltatol)]<-0
  
  ### $\tfrac{1}{2}\sum_{j=1}^{p}\hat{\mu}_{kj}^2$ constant term
  consts<-rep(0,K)
  for(k in 1:K) consts[k]<-0.5*sum(MuIter[k,]^2/sigmasqs)
  ### $\hat{\mu}_{kj}^2/\sigma_p^2$ term
  lins<-vector(mode="list",length=K)
  for(k in 1:K) lins[[k]]<-MuIter[k,]/sigmasqs
  
  ### return calculated information as list object
  return(list(classes=classnames,consts=consts,lins=lins,prior=priors
              ,meanadj=transMeans,MuIter=MuIter,InitMu=MuMat))
}
