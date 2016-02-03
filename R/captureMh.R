#'This function simulates the capture process with a 2-component finite mixture model for individual heterogeneity.
#'@param p a vector of size 3 containing the two capture probabilities and the mixture proportion  p=c(p1,p2,mix)
#'@param N an integer that is the true population size
#'@param occ an integer that is the number of capture occasions to simulate
#'@param traptype a character string specifying the trap type.  "single" 
#'allows individuals to be caught in at most 1 trap per occasion and "multi" allows
#'individuals to be caught in multiple traps per occasion.
#'@param K an integer that is the number of traps to simulate if traptype="multi"
#'@return W a capture history of dimension N x occ if traptype="single" or 
#'N x occ x K if traptype="multi"
#'@export
#'
captureMh=function(p,N,occ,traptype,K=NULL){
  if(length(p)!=3){stop("p must be of dimension 3")}
  #get heterogeneity
  pmix=rbinom(N,1,p[3])
  phet=ifelse(pmix==1,p[1],p[2])
  if(traptype=="single"){
    #Simulate capture process with single trap captures per occasion and hair deposition
    W=matrix(rbinom(N*occ,1,phet),ncol=occ,nrow=N)
  }else if(traptype=="multi"){
    if(is.null(K)){
      stop("if traptype=multi, you must enter the number of traps, K")
    }
    #Simulate capture process with multiple trap captures per occasion and hair deposition
    ptrap=1-(1-phet)^(1/K)
    W=array(0,dim=c(N,occ,K)) #latent capture history
      for(j in 1:occ){
        W[,j,]=matrix(rbinom(N*K,1,ptrap),ncol=K,nrow=N)
      }
  }
  return(W)
}

