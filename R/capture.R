#'This function simulates the capture process
#'@param p a vector of size 2 containing the capture and recapture probabilities
#'@param N an integer that is the true population size
#'@param occ an integer that is the number of capture occasions to simulate
#'@param traptype a character string specifying the trap type.  "single" 
#'allows individuals to be caught in at most 1 trap per occasion and "multi" allows
#'individuals to be caught in multiple traps per occasion.
#'@param btype a character string specifying the type of behavioral response to capture
#'Current options are "global" for a global trap response and "trap" for a trap-specific
#'response.  "trap" currently only works with traptype=multi.
#'@param K an integer that is the number of traps to simulate if traptype="multi"
#'@return W a capture history of dimension N x occ if traptype="single" or 
#'N x occ x K if traptype="multi"
#'@export
#'
capture=function(p,N,occ,traptype,btype,K=NULL){
  if(traptype=="single"){
    if(btype=="trap"){
      stop("Not implementing btype=trap if traptype=single")
    }
    state=rep(1,N) #Animal states.  1=not previously captured 2=previously captured
    #Simulate capture process with single trap captures per occasion and hair deposition
    W=matrix(0,nrow=N,ncol=occ) #latent capture history
    for(i in 1:occ){
      W[,i]=rbinom(N,1,p[state])
      state[W[,i]==1 & state==1]=2
    }
  }else if(traptype=="multi"){
    if(is.null(K)){
      stop("if traptype=multi, you must enter the number of traps, K")
    }
    #Simulate capture process with multiple trap captures per occasion and hair deposition
    ptrap=1-(1-p)^(1/K)
    W=array(0,dim=c(N,occ,K)) #latent capture history
    if(btype=="global"){
      state=rep(1,N) #Animal states.  1=not previously captured 2=previously captured
      for(j in 1:occ){
        W[,j,]=matrix(rbinom(N*K,1,ptrap[state]),ncol=K,nrow=N)
        caught=W[,j,]==1
        caught2=(rowSums(caught)>0)*1
        state[caught2==1 & state==1]=2
      }
    } else if (btype=="trap"){
      state=matrix(1,nrow=N,ncol=K) #Animal states.  1=not previously captured 2=previously captured
      for(j in 1:occ){
        W[,j,]=matrix(rbinom(N*K,1,ptrap[state]),ncol=K,nrow=N)
        caught=W[,j,]==1
        state[caught==1 & state==1]=2
      }
    }
  }
  return(W)
}