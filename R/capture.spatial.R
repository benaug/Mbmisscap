#'This function simulates a spatial capture process
#'@param N an integer that is the true population size
#'@param occ an integer that is the number of capture occasions to simulate
#'@param lambda
#'@param sigma
#'@param locs
#'@param buff
#'@param btype a character string specifying the type of behavioral response to capture
#'Current options are "global" for a global trap response and "trap" for a trap-specific
#'response.  "trap" currently only works with traptype=multi.
#'@return W a capture history of dimension N x occ x K ...
#'@export
#'
capture.spatial=function(N,occ,lambda,sigma,locs,buff,btype=NULL){
  #Get buffer limits for simulating activity centers
  lowerx=min(locs[,1])-buff
  upperx=max(locs[,1])+buff
  lowery=min(locs[,2])-buff
  uppery=max(locs[,2])+buff
  #Simulate N home range centers
  ACx=runif(N, lowerx,upperx)
  ACy=runif(N, lowery,uppery)
  #Make the Dij matrix of distances from home center to trap locations
  K=nrow(locs)
  D=matrix(0, nrow=N, ncol=K)
  for(i in 1:N){
    D[i,] = spDistsN1(as.matrix(locs), c(ACx[i], ACy[i]))
  }
  #k=exp(-D^2/sigma) Royle uses this
  #Exposure to traps g
  g=exp(-D^2/(2*sigma^2))
  pmean=matrix(NA, nrow=N, ncol=K)
  cmean=matrix(NA, nrow=N, ncol=K)
  for(i in 1:N){
    for(k in 1:K){
      pmean[i,k] = 1 -exp(-lambda[1]*g[i,k]) 
      cmean[i,k] = 1 -exp(-lambda[2]*g[i,k])       
    }
  }
  #trap specific Animal states.  1=not previously captured 2=previously captured
  qmean=abind(pmean,cmean,along=3)
  W=array(0,dim=c(N,occ,K)) #latent capture history
  if(btype=="trap"){
    state=matrix(1,nrow=N,ncol=K) 
    for(i in 1:N){
      for(j in 1:occ){
        for(k in 1:K){
          #check this
          W[i,j,k]=rbinom(1,1,qmean[i,k,state[i,k]])
          if(state[i,k]==1&W[i,j,k]==1){
            state[i,k]=2
          }
        }
      }
    }  
  } else {  #add option for "none" instead of reusing this?
    state=rep(1,N) 
    for(i in 1:N){
      for(j in 1:occ){
        for(k in 1:K){
          #check this
          W[i,j,k]=rbinom(1,1,qmean[i,k,state[i]])
          if(state[i]==1&W[i,j,k]==1){
            state[i]=2
          }
        }
      }
    }  
  }
  p_i=apply(pmean,1,FUN=function(x){1-prod(1-x)})
  c_i=apply(cmean,1,FUN=function(x){1-prod(1-x)})
  #Record activity centers and order them to match capture history
  ACs=data.frame(x=ACx,y=ACy)
  ACs$cap=1*(apply(W,1,sum)>0)
  return(list(W=W,p_i=p_i,c_i=c_i,ACs=ACs))
}