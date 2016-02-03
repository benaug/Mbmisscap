#'This function simulates a hair snare mark recapture experiment with hair sample subsampling and 
#'failed DNA identification.
#'@param N an integer that is the true population size
#'@param p a vector of size 2 containing the capture and recapture probabilities...fix this
#'@param occ an integer that is the number of capture occasions to simulate
#'@param K an integer that is the number of traps to simulate
#'@param traptype a character string specifying the trap type.  "single" 
#'allows individuals to be caught in at most 1 trap per occasion and "multi" allows
#'individuals to be caught in multiple traps per occasion.
#'@param lambda_h a positive value specifying the parameter for the zero-truncated Poisson
#'hair deposition process (S_ijk|W_ijk=1)
#'@param lambda_c a positive value specifying the parameter for the zero-truncated Poisson
#'cluster deposition process if cluster=TRUE
#'@param kappa_h an integer specifying the number of hair samples to retain per occasion-trap if cluster=FALSE
#'or per occasion-trap-cluster if cluster=TRUE.  Used in sub-hair and sub-cluster methods
#'@param kappa_t an integer specifying the number of traps per occasion at which to retain kappa_h hair samples.
#'If no value is specified no trap-level subsampling occurs.  Only works when cluster=FALSE.  Used in the
#'sub-trap method
#'@param delta a numeric value between 0 and 1 specifying the probability a hair sample will be retained
#'in the subsample.  This is used if subsampling is done by pooling samples on each occasion and taking
#'a simple random sample.
#'@param cluster a logical indicating whether or not to simulate from the cluster model (multiple clusters per
#'individual-occasion-trap)
#'@param kappa_c an integer specifying the number of clusters per trap-occasion at which to retain kappa_h hair samples.
#'Only works when cluster=TRUE.  Sub-cluster method.
#'@param alpha a numeric value between 0 and 1 specifying the probability that a hair sample produces
#'an individual identification
#'param hettype a character string specifying the type of individual heterogeneity.  Either "logitnormal" or "finitemixture".
#'@return a list with the simulated quantities and statistics (expand description later)
#'@include captureMh2
#'@include deposit
#'@include subsample
#'@include DNAamp
#'@include calcSS
#'@include calcSUR
#'@export
simMhmisscap=function(N,p,occ,K,traptype,lambda_h,lambda_c=NULL,
                  delta=NULL,kappa_h=NULL,kappa_t=NULL,cluster=FALSE,kappa_c=NULL,alpha=NULL,hettype="logitnormal"){
  ##check for input errors  
  if(is.null(kappa_h)&is.null(delta)){
    stop("Must enter either kappa_h or delta")
  }
  if(!is.null(kappa_h)&!is.null(delta)){
    stop("Cannot enter both kappa_h and delta")
  }
  if(cluster==TRUE&is.null(lambda_c)){
    stop("if cluster=TRUE, must specify lambda_c")
  }
  if(hettype=="logitnormal"&length(p)==3){
    stop("If hettype==logitnormal, p must be of dimension 2")
  }
  if(hettype=="finitemixture"&length(p)==2){
    stop("If hettype==finitemixture, p must be of dimension 3")
  }
  ##Capture process
  if(hettype=="logitnormal"){
    W=captureMh2(p,N,occ,traptype,K)    
  } else {
    W=captureMh(p,N,occ,traptype,K)
  }
  
  ##Hair deposition process
  S=deposit(W,lambda_h,traptype,K,cluster,lambda_c)
  ##Subsampling process
  Sdot=apply(S,2:3,sum)
  U=subsample(S,Sdot,delta=delta,kappa_h=kappa_h,kappa_t=kappa_t,cluster=cluster,kappa_c=kappa_c)
  ##DNA Identification process
  if(!is.null(alpha)){
    Udot=apply(U,2:3,sum)
    R=DNAID(alpha,U)
  }else{
    Udot=Sdot
    R=U
  }

  #######build observed capture history######
  if(cluster==TRUE){ #just collapse for now to make compatible
    R=apply(R,1:3,sum)
    S=apply(S,1:3,sum)
    U=apply(U,1:3,sum)
    W=(apply(W,1:3,sum)>0)*1
  }
  Wobsfull=(R>0)*1
  empty=apply(Wobsfull,1,function(x){all(x==0)})
  if(any(empty==TRUE)){ 
    Wobs=Wobsfull[-which(empty==T),,]
    S=S[-which(empty==TRUE),,]
    U=U[-which(empty==TRUE),,]
    R=R[-which(empty==TRUE),,]
  } else {
    Wobs=Wobsfull
  }
  n=dim(Wobs)[1]
  
  ####Calculate capture history stats#####
  #Collapse capture histories
  W2D=(apply(W,1:2,sum)>0)*1
  Wobs2D=(apply(Wobs,1:2,sum)>0)*1
  Wobsfull2D=(apply(Wobsfull,1:2,sum)>0)*1
  Wstats=calcSS(W2D,Wobs2D,Wobsfull2D)

  ####Calculate hair history stats#####
  S2D=apply(S,1:2,sum)
  U2D=apply(U,1:2,sum)
  R2D=apply(R,1:2,sum)

  HairStats=calcSUR(S2D,U2D,R2D)
  return(list(n=n,
              occ=occ,
              K=K,
              Wobs=Wobs,
              W=W,
              W2D=W2D,
              Wobs2D=Wobs2D,
              S=S,
              U=U,
              R=R,
              Sdot=Sdot,
              Udot=Udot,
              Wstats=Wstats,
              HairStats=HairStats))
}