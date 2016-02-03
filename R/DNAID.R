#'This function simulates the DNA identification process
#'@param U a "subsampled hair capture history" of dimension N x occ x K if cluster=FALSE or 
#'N x occ x K x Lmax if cluster=TRUE
#'@param alpha a numeric value between 0 and 1 specifying the probability that a hair sample produces
#'an individual identification
#'@return R an "observed hair capture history" containing the number of hair samples deposited by individual i on occasion
#'j in trap k and if for cluster=TRUE, at cluster l, that are retained in the subsample and produce an individual
#'identification.  R is of dimension N x occ x K if cluster=FALSE and dimension N x occ x K x Lmax if cluster=TRUE,
#' where Lmax is the maximum observed number of clusters deposited across trap-occasions
#'@export
#'
DNAID=function(alpha,U){
    R=array(rbinom(prod(dim(U)),U,alpha),dim=dim(U))
    return(R)
}