####### FUNCTION: PFldapredict()
### estimate probabilities and class of new inputs

####### input:
### ldaobj: an object created by createPFldaobj()
### Xnew: a n x p matrix, of n obs and p variables. 
###         May also be a single numeric vector (n=1) of length p
### priors: a vector of length K (#classes) with elements in (0,1)

PFldapredict<-function(ldaobj,Xnew,priors=NULL)
{
  if(is.vector(Xnew,mode="numeric")){
    Xnew<-matrix(Xnew,nrow=1)
  }else if(!is.matrix(Xnew)){
    cat("Xnew must be numeric and either a vector or matrix \n ")
  }

  Nnew<-nrow(Xnew)
  ### centre each feature at 0, as done on the training data
  Xnew<-Xnew-matrix(rep(ldaobj$meanadj,Nnew),nrow=Nnew,byrow=TRUE)
  Xclasses<-ldaobj$classes
  K<-length(Xclasses)
  if(is.null(priors)) priors<-ldaobj$prior
  
  outposteriors<-matrix(0,nrow=Nnew,ncol=K)
  colnames(outposteriors)<-Xclasses
  
  ### discrim function: $\delta_k(x_{new})$
  for(i in 1:Nnew){
    for(k in 1:K){
      outposteriors[i,k]<-log(priors[k]) - ldaobj$consts[k] 
                                 + sum(Xnew[i,]*ldaobj$lins[[k]])
    }
  }
  
  ### predicted class is arg max
  ### $p(G_K|x_{new})=\tfrac{P\left(x_{new}|G_k\right)P\left(G_k\right)}{\sum_{i=1}^{K}P\left(x_{new}|G_i\right)P\left(G_i\right)}$
  predclasses<-rep(0,Nnew)
  for(i in 1:Nnew){
    outposteriors[i,]<-exp(outposteriors[i,])
    predclasses[i]<-which.max(outposteriors[i,])
    outposteriors[i,]<-outposteriors[i,]/sum(outposteriors[i,])
  }
  outpred<-factor(Xclasses[predclasses])
  
  return(list(pred=outpred,posteriors=outposteriors))
}