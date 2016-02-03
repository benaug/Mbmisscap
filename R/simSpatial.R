#'This function simulates a spatial hair snare mark recapture experiment with hair sample subsampling and 
#'failed DNA identification.
#'@param N an integer that is the true population size
#'@param occ an integer that is the number of capture occasions to simulate
#'@param lambda
#'@param sigma
#'@param locs
#'@param buff
#'@param btype a character string specifying the type of behavioral response to capture
#'Current options are "global" for a global trap response and "trap" for a trap-specific
#'response.  "trap" currently only works with traptype=multi.
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
#'@return a list with the simulated quantities and statistics (expand description later)
#'@include capture
#'@include deposit
#'@include subsample
#'@include DNAamp
#'@include calcSS
#'@include calcSUR
#'@export
simSpatial=function(N,occ,lambda,sigma,locs,buff,btype,lambda_h,lambda_c=NULL,
                          delta=NULL,kappa_h=NULL,kappa_t=NULL,cluster=FALSE,kappa_c=NULL,alpha=NULL){
  ##check for input errors  
  if(is.null(kappa_h)&is.null(delta)){
    warning("Must enter either kappa_h or delta")
  }
  if(is.null(delta)&is.null(kappa_h)){
    stop("Must enter either delta or kappa_h")
  }
  if(!is.null(kappa_h)&!is.null(delta)){
    warning("Cannot enter both kappa_h and delta")
  }
  require(sp)
  require(abind)
  require(secr)
  ##Capture process
  cap=capture.spatial(N,occ,lambda,sigma,locs,buff,btype)
  W=cap$W
  
  ##Hair deposition process
  K=nrow(locs)
  S=deposit(W,lambda_h,K=K,cluster=cluster,lambda_c=lambda_c)
  
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
  
  #Update ACs after subsampling
  cap$ACs$cap2=1*(apply(Wobsfull,1,sum)>0)
  
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
              AC=cap$AC,
              Wstats=Wstats,
              HairStats=HairStats,
              pstats=Wstats$pstats
             ))
}