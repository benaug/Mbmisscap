#'This function simulates the subsampling process
#'@param S a "hair capture history" of dimension N x occ x K if cluster=FALSE or 
#'N x occ x K x Lmax if cluster=TRUE
#'@param Sdot a matrix of constraints on S obtained by summing S across the second and third dimensions.
#'apply(S,2:3,sum)  Sdot will be 2 dimensional if cluster=FALSE and 3 dimensional if cluster=TRUE
#'@param delta a numeric value between 0 and 1 specifying the probability a hair sample will be retained
#'in the subsample.  This is used if subsampling is done by pooling samples on each occasion and taking
#'a simple random sample.
#'@param kappa_h an integer specifying the number of hair samples to retain per occasion-trap if cluster=FALSE
#'or per occasion-trap-cluster if cluster=TRUE.  Used in sub-hair and sub-cluster methods
#'@param kappa_t an integer specifying the number of traps per occasion at which to retain kappa_h hair samples.
#'If no value is specified no trap-level subsampling occurs.  Only works when cluster=FALSE.  Used in the
#'sub-trap method
#'@param cluster a logical indicating whether or not to simulate from the cluster model (multiple clusters per
#'individual-occasion-trap)
#'@param kappa_c an integer specifying the number of clusters per trap-occasion at which to retain kappa_h hair samples.
#'Only works when cluster=TRUE.  Sub-cluster method.
#'@return U a "subsampled hair capture history" containing the number of hair samples deposited by individual i on occasion
#'j in trap k and if for cluster=TRUE, at cluster l, that are retained in the subsample.  U is of dimension N x occ x K if cluster=FALSE and dimension
#'N x occ x K x Lmax if cluster=TRUE, where Lmax is the maximum observed number of clusters deposited across
#'trap-occasions
#'@export
#'
subsample=function(S,Sdot,delta=NULL,kappa_h=NULL,kappa_t=NULL,cluster=FALSE,kappa_c=NULL){
  #Some user error checks
  if(!is.null(delta)&!is.null(kappa_h)){
    stop("Must specify either delta or kappa_h")
  }
  if(is.null(delta)&is.null(kappa_h)){
    stop("Must specify either delta or kappa_h")
  }
  if(!is.null(kappa_t)&is.null(kappa_h)){
    stop("if you specify kappa_t, you must specify kappa_h")
  }
  if(!is.null(kappa_c)&is.null(kappa_h)){
    stop("if you specify kappa_c, you must specify kappa_h")
  }
  
  N=dim(S)[1]
  occ=dim(S)[2]
  K=dim(S)[3]
  if(cluster==FALSE){
    if(is.null(delta)){  #Sub methods.  No SRS.
      if(is.null(kappa_t)){   #Sub-hair method
        U=S #Set U=S to fix U where S[,j,k]<kappa_h.  Subsampling below replaces values in U where subsampling is done
        idx=which(Sdot>kappa_h,arr.ind=TRUE) #which ind-occ-trap have >kappa_h hair samples?
        for(k in 1:nrow(idx)){
          U[,idx[k,1],idx[k,2]]=rMVhyper(n=kappa_h,K=S[,idx[k,1],idx[k,2]])
        }
      #Sub-trap method
      }else{
        U=S*0
        for(j in 1:occ){
          #Select subsample of traps.  
          usetrap=which(Sdot[j,]>0)#Start with those that have at least 1 hair sample          
          if(length(usetrap)>kappa_t){
            usetrap=sample(usetrap,kappa_t) #subsample traps w/o replacement
          }
          #set U=S for traps w/o subsampling
          U[,j,usetrap]=S[,j,usetrap]
          for(k in usetrap){
            if(Sdot[j,k]>kappa_h){
              U[,j,k]=rMVhyper(n=kappa_h,K=S[,j,k])
            }
          }
        }
      }
      #SRS
    }else{
      U=array(rbinom(N*occ*K,S,delta),dim=dim(S))
    }
  }else{ #Sub-cluster method
    if(is.null(kappa_c)){  #no cluster subsampling
    U=S #Set U=S and fix U where S[,j,k]<kappa and subsampling must be done
    idx=which(S>kappa_h,arr.ind=TRUE) #which ind-occ-trap-clusters have >kappa hair samples?
    #MVhyper works with or without contamination.  Without, it defaults to hypergeometric.
    for(l in 1:nrow(idx)){
      U[,idx[l,2],idx[l,3],idx[l,4]]=rMVhyper(n=kappa_h,K=S[,idx[l,2],idx[l,3],idx[l,4]])
    }
    } else { #Subsample clusters
      U=S*0
      Sdot2=apply(S,2:4,sum)
      for(j in 1:occ){
        trapidx=which(Sdot[j,]>1) #which traps on each occasion have at least one cluster?
        for(k in trapidx){
          clustidx=which(Sdot2[j,k,]>0) #id clusters at this trap-occasion
          if(length(clustidx)>kappa_c){ #do we need to subsample clusters?
            clustkeep=sample(clustidx,kappa_c)
          } else{
            clustkeep=clustidx
          }
          for(l in clustkeep){
            U[,j,k,l]=rMVhyper(n=kappa_h,K=S[,j,k,l])
          }
        }
      }
    }
  }
  return(U=U)
}