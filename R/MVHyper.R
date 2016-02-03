#'PDF for the multivariate hypergeometric distribution
#'@param k adf
#'@param K adf
#'@export
#'

dMVhyper=function(k,K){
  n=sum(k)
  N=sum(K)
  p=prod(choose(K,k))/choose(N,n)
  return(p)
}
#'Random generator for the multivariate hypergeometric distribution
#'@param K adf
#'@param n ads
#'@export
#'
rMVhyper=function(n,K){
  vect=rep(1:length(K),times=K)
  if(sum(K)>0&length(vect)>1){
    N=sum(K)
    samp=sample(vect,n)
    a=1:length(K)
    k=sapply(a,function(x){sum(a[x]==samp)})
  }else{
    k=K
  }
  return(k)
}

