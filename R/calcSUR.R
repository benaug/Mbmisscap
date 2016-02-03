#'This function calculates the failure criterion for Mb (see Otis 1978) from the observed capture history Wobs2D
#'and the sufficient statistics for Mb from both the true and observed capture histories.  (add description of SS later)
#'@param S the true hair capture history of dimension N x occ
#'@param U the hair capture history after subsampling of dimension N x occ
#'@param R the hair capture history after subsampling and DNA Identificationof dimension N x occ
#'@return a list containing the stuff
#'@export
#'
calcSUR=function(S2D,U2D,R2D){
  Sdotj=colSums(S2D)
  Udotj=colSums(U2D)
  Rdotj=colSums(R2D)
  return(list(Sdotj=Sdotj,
              Udotj=Udotj,
              Rdotj=Rdotj))
}