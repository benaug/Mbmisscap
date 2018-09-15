#'This function simulates the hair deposition process
#'@param W a capture history of dimension N x occ if traptype="single" or
#'N x occ x K if traptype="multi"
#'@param lambda_h a positive value specifying the parameter for the zero-truncated Poisson
#'hair deposition process (S_ijk|W_ijk=1)
#'@param traptype a character string specifying the capture process.  "single"
#'allows individuals to be caught in at most 1 trap per occasion and "multi" allows
#'individuals to be caught in multiple traps per occasion.  If traptype="single", a trap id is randomly
#'assigned to each capture
#'@param K an integer that is the number of traps to simulate
#'@param cluster a logical indicating whether or not to simulate from the cluster model (multiple clusters per
#'individual-occasion-trap)
#'@param lambda_c a positive value specifying the parameter for the zero-truncated Poisson
#'cluster deposition process for the cluster model
#'@return S a "hair capture history" containing the number of hair samples deposited by individual i on occasion
#'j in trap k and if for cluster=TRUE, at cluster l.  S is of dimension N x occ x K if cluster=FALSE and dimension
#'N x occ x K x Lmax if cluster=TRUE, where Lmax is the maximum observed number of clusters deposited across
#'trap-occasions
#'@import VGAM
#'@export
#'
deposit=function(W,lambda_h,traptype="multi",K,cluster=FALSE,lambda_c=NULL){
  N=dim(W)[1]
  occ=dim(W)[2]
  if(cluster==FALSE){
    S=array(0,dim=c(N,occ,K)) #deposited hair samples
    if(traptype=="single"){
      #Simulate hair depostion with 2D W
      for(i in 1:occ){
        for(j in 1:N){
          if(W[j,i]>0){
            #Randomly assign a trap at random if caught
            S[j,i,sample(1:K,1)]=VGAM::rzapois(1,lambda_h,pobs0=0)
          }
        }
      }
    }else if(traptype=="multi"){
      #simulate hair deposition with 3D W
      for(i in 1:occ){
        hair=matrix(VGAM::rzapois(N*K,lambda_h,pobs0=0),ncol=K,nrow=N)
        caught=W[,i,]==1
        S[,i,][caught]=hair[caught]
      }
    }
  }else{
    #cluster deposition model, multiple traps only at the moment
    C=W*0
    C[W==1]=VGAM::rzapois(sum(W),lambda_c,pobs0=0) #clusters per bear-trap-occasion
    #Lmax=max(C)
    Cdot=apply(C,2:3,sum)
    Lmax=max(Cdot)
    #Build D matrix
    D=array(0,dim=c(N,occ,K,Lmax)) #cluster indicator matrix
    idx=which(Cdot>0,arr.ind=TRUE) #identify occasion-traps with clusters to assign
    for(l in 1:nrow(idx)){
      Cidotdot=C[,idx[l,1],idx[l,2]]
      inds=which(Cidotdot>0)
      Cinddotdot=Cidotdot[inds]
      clus=1
      for(i in 1:length(inds)){
        D[inds[i],idx[l,1],idx[l,2],clus:(clus+Cinddotdot[i]-1)]=rep(1,Cinddotdot[i])
        clus=clus+Cinddotdot[i]
      }
    }

    #Build D matrix ind-specific cluster IDs - old
#     D=array(0,dim=c(N,occ,K,Lmax)) #cluster indicator matrix
#     idx=which(C>0,arr.ind=TRUE) #which ind-occ-traps have least one cluster?
#     for(l in 1:nrow(idx)){
#       lclus=C[idx[l,1],idx[l,2],idx[l,3]]
#       D[idx[l,1],idx[l,2],idx[l,3],1:lclus]=1
#     }


    S=array(0,dim=c(N,occ,K,Lmax)) #deposited hair samples at each cluster
    S[D==1]=VGAM::rzapois(sum(D),lambda_h,pobs0=0) #clusters per bear-trap-occasion
  }
  return(S=S)
}
