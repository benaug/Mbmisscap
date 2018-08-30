#'This function calculates the failure criterion for Mb (see Otis 1978) from the observed capture history Wobs2D
#'and the sufficient statistics for Mb from both the true and observed capture histories.  (add description of SS later)
#'@param W the true capture history of dimension N x occ
#'@param Wobs2D the observed capture history (after subsampling and DNA amplification) of dimension N x occ
#'@return a list containing the sufficient statistics and failure criterion
#'@export
#'
calcSS=function(W2D,Wobs2D,Wobsfull2D){
  occ=ncol(W2D)

  #######Calculate missing first caps
  firstTrue=rep(0,nrow(W2D))
  for(i in 1:nrow(W2D)){
    firstTrue[i]=match(1,W2D[i,])
  }
  firstcapsTrue=as.numeric(table(firstTrue))
  firstObsfull=rep(0,nrow(Wobsfull2D))
  for(i in 1:nrow(Wobsfull2D)){
    firstObsfull[i]=match(1,Wobsfull2D[i,])
  }
  miss=table(firstTrue,(firstTrue<firstObsfull|is.na(firstObsfull))*1)  #1 if first cap missing 0 ow
  miss1caps=miss[,2]
  diff=occ-length(miss1caps)
  if(diff>0){
    miss1caps=c(miss1caps,rep(0,diff))
  }
  diff=occ-length(firstcapsTrue)
  if(diff>0){
    firstcapsTrue=c(firstcapsTrue,rep(0,diff))
  }
  ##Calculate p stats
  firstTruefact=factor(firstTrue)
  firstObsfullfact=factor(firstObsfull)
  levels(firstTruefact)=1:occ
  levels(firstObsfullfact)=1:occ
  tab=table(firstTruefact,firstObsfullfact)
  totalps=colSums(tab)
  tab2=tab
  diag(tab2)=0
  nfalsecaps=colSums(tab2)
  pinpoptrue=c(N,(N-cumsum(table(firstTrue)))[1:(occ-1)])
  pinpopobs=c(N,(N-cumsum(table(firstObsfull)))[1:(occ-1)])
  pfalseinpop=pinpopobs-pinpoptrue
  pfalse=nfalsecaps/pfalseinpop #capture probability of falsely identified first captures
  pfalse[1]=0
  ptrue=diag(tab)/pinpoptrue # true capture probability of correctly identified first captures
  ptrueobs=diag(tab)/pinpopobs # estimated capture probability of correctly identified first captures
  perfalse=nfalsecaps/totalps # percent of identified first captures that are false
  pobs=pfalse*perfalse+ptrueobs*(1-perfalse)

  pstats=rbind(ptrue,pobs,ptrueobs,pfalse,perfalse)
  rownames(pstats)=c("p_true","p_obs","p_obs_true","p_obs_false","percent_false")
  ####Calculate True Sufficient Statistics#####
  W2D=W2D[rowSums(W2D)>0,] #remove true all zero capture histories and recalculate some things from above

  #m_js
  firstTrue=rep(0,nrow(W2D))
  for(i in 1:nrow(W2D)){
    firstTrue[i]=match(1,W2D[i,])
  }
  mj_true=rep(0,occ)
  for(i in 2:occ){
    mj_true[i]=length(which(firstTrue[which((W2D[,i]==1))]<i))
  }

  #M_js
  Mj_true=c(0,cumsum(table(firstTrue))[-ncol(W2D)])
  diff=ncol(W2D)-length(Mj_true)
  if(diff>0){
    Mj_true=c(Mj_true,rep(Mj_true[length(Mj_true)],diff))
  }

  #mdot, Mdot, and Mt+1
  mdot_true=sum(mj_true)
  Mdot_true=sum(Mj_true)
  Mt1_true=nrow(W2D)

  ####Calculate Observed Sufficient Statistics######
  #m_js
  firstObs=rep(0,nrow(Wobs2D))
  for(i in 1:nrow(Wobs2D)){
    firstObs[i]=match(1,Wobs2D[i,])
  }
  firstcapsObs=as.numeric(table(firstObs))##
  mj_obs=rep(0,occ)
  for(i in 2:occ){
    mj_obs[i]=length(which(firstObs[which((Wobs2D[,i]==1))]<i))
  }

  #M_js
  Mj_obs=c(0,cumsum(table(firstObs))[-ncol(Wobs2D)])
  diff=ncol(W2D)-length(Mj_obs)
  if(diff>0){
    Mj_obs=c(Mj_true,rep(Mj_obs[length(Mj_obs)],diff))
  }

  #mdot, Mdot and Mt+1
  mdot_obs=sum(mj_obs)
  Mdot_obs=sum(Mj_obs)
  Mt1_obs=nrow(Wobs2D)

  #Observed failure criterion
  fail=sum((occ+1-2*1:occ)*(colSums(Wobs2D)-mj_obs))
  return(list(mj_true=mj_true,
              Mj_true=Mj_true,
              mdot_true=mdot_true,
              Mdot_true=Mdot_true,
              Mt1_true=Mt1_true,
              mj_obs=mj_obs,
              Mj_obs=Mj_obs,
              mdot_obs=mdot_obs,
              Mdot_obs=Mdot_obs,
              Mt1_obs=Mt1_obs,
              firstcapsTrue=firstcapsTrue,
              miss1caps=miss1caps,
              pstats=pstats,
              fail=fail))
}
