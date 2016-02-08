#'This function asseses some frequentist properties when naively fitting Mb to capture histories that are missing data.
#'@param N an integer that is the true population size
#'@param p a vector of size 2 containing the capture and recapture probabilities
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
#'@param hettype ...
#'@param sims the number of simulated data sets to generate and fit Mb to
#'@param cores the number of cores to do the simulation on.  A typical for loop is used if cores=1 and
#'the foreach package is used if cores>1
#'@return a list with the simulation results (expand later)
#'@export

assessMh=function(N,p,occ,K,traptype,lambda_h,lambda_c=NULL,delta=NULL,kappa_h=NULL,
                kappa_t=NULL,cluster=FALSE,kappa_c=NULL,alpha=NULL,hettype="logitnormal",sims=100,cores=1){

  if(hettype=="logitnormal"){
    if(length(p)==3){
      stop("If hettype==logitnormal, p must be of dimension 2")
    }
    #require(Rcapture)
  }
  if(hettype=="finitemixture"){
    if(length(p)==2){
      stop("If hettype==finitemixture, p must be of dimension 3")
    }
  }
  if(cores==1){
    #preallocate storage structures
    storeN=matrix(NA,nrow=sims,ncol=3)
    storeCapsTrue=matrix(NA,nrow=sims,ncol=occ)
    storeCapsObs=matrix(NA,nrow=sims,ncol=occ)
    storeS=matrix(NA,nrow=sims,ncol=occ)
    storeU=matrix(NA,nrow=sims,ncol=occ)
    storeR=matrix(NA,nrow=sims,ncol=occ)
    for(i in 1:sims){
      if(hettype=="logitnormal"){
        data=simMhmisscap(N,p,occ,K,traptype,lambda_h,lambda_c,delta,kappa_h,
                          kappa_t,cluster,kappa_c,alpha,hettype="logitnormal")
        mod=Rcapture::closedpCI.0(data$Wobs2D, m = "Mh", h = "Normal")
        storeN[i,]=mod$results[c(1,3,4)]
      } else{
        data=simMhmisscap(N,p,occ,K,traptype,lambda_h,lambda_c,delta,kappa_h,
                          kappa_t,cluster,kappa_c,alpha,hettype="finitemixture")
        Mdata=as.data.frame(apply(data$Wobs2D,1,paste,collapse=""))
        colnames(Mdata)="ch"
        Mdata$ch=as.character(Mdata$ch)
        out=RMark::mark(Mdata,model="FullHet",model.parameters=list(p=list(formula=~mixture,share=TRUE)),output=FALSE,delete=TRUE) #Mh
        storeN[i,]=unlist(out$results$derived[c(1,3,4)])
      }
      storeCapsTrue[i,]=colSums(data$W2D)
      storeCapsObs[i,]=colSums(data$Wobs2D)
      storeS[i,]=data$HairStats$Sdotj
      storeU[i,]=data$HairStats$Udotj
      storeR[i,]=data$HairStats$Rdotj
    }
    den=sims
  }else{
    #load multicore packages
    require(foreach)
    require(snow)
    require(doSNOW)
    cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp)
    clusterExport(cl.tmp, list=ls(environment()), envir = environment())
    if(hettype=="logitnormal"){
      require(Rcapture)
      out2=foreach(i=1:sims,.packages=c("VGAM","Rcapture","Mbmisscap"),.errorhandling = "remove") %dopar% {
        data=simMhmisscap(N,p,occ,K,traptype,lambda_h,lambda_c,delta,kappa_h,
                          kappa_t,cluster,kappa_c,alpha,hettype="logitnormal")
        mod=closedpCI.0(data$Wobs2D, m = "Mh", h = "Normal")

        Nhat=mod$results[c(1,3,4)]
        return(list(Nhat,
                    colSums(data$W2D),
                    colSums(data$Wobs2D),
                    data$HairStats$Sdotj,
                    data$HairStats$Udotj,
                    data$HairStats$Rdotj))
      }
    } else {
      require(RMark)
      out2=foreach(i=1:sims,.packages=c("VGAM","RMark"),.errorhandling = "remove") %dopar% {
        data=simMhmisscap(N,p,occ,K,traptype,lambda_h,lambda_c,delta,kappa_h,
                          kappa_t,cluster,kappa_c,alpha,hettype="finitemixture")
        Mdata=as.data.frame(apply(data$Wobs2D,1,paste,collapse=""))
        colnames(Mdata)="ch"
        Mdata$ch=as.character(Mdata$ch)
        out=mark(Mdata,model="FullHet",model.parameters=list(p=list(formula=~mixture,share=TRUE)),output=FALSE,delete=TRUE) #Mh
        f0=get.real(out,"f0",se=TRUE)[3:4]
        c=exp(1.96*sqrt((log(1+f0[2]^2/(f0[1]^2)))))
        Nhat=unlist(c(data$n+f0[1],data$n+f0[1]/c,data$n+f0[1]*c))

        return(list(Nhat,
                    colSums(data$W2D),
                    colSums(data$Wobs2D),
                    data$HairStats$Sdotj,
                    data$HairStats$Udotj,
                    data$HairStats$Rdotj))
      }
      cleanup(ask=FALSE,lx=ls(.GlobalEnv))
    }
    stopCluster(cl.tmp)
    den=length(out2) #successful fits RMark didn't loose track of
    #preallocate storage structures
    storeN=matrix(NA,nrow=den,ncol=3)
    storeFail=rep(NA,den)
    storeCapsTrue=matrix(NA,nrow=den,ncol=occ)
    storeCapsObs=matrix(NA,nrow=den,ncol=occ)
    storeS=matrix(NA,nrow=den,ncol=occ)
    storeU=matrix(NA,nrow=den,ncol=occ)
    storeR=matrix(NA,nrow=den,ncol=occ)
    for(i in 1:den){
      storeN[i,]=out2[[i]][[1]]
      storeCapsTrue[i,]=out2[[i]][[2]]
      storeCapsObs[i,]=out2[[i]][[3]]
      storeS[i,]=out2[[i]][[4]]
      storeU[i,]=out2[[i]][[5]]
      storeR[i,]=out2[[i]][[6]]
    }
  }



  #Cacluate N stats
  Nmean=mean(storeN[,1],na.rm=TRUE)
  dir=Nmean[1]-N
  Nbias=ifelse(dir>0,paste("+",round(100*(Nmean[1]/N-1),2),"%"),paste("-",round(100*(1-Nmean[1]/N),2),"%"))
  Ncover=round(sum(storeN[,2]<N&storeN[,3]>N,na.rm=TRUE)/den,2)
  Nwidth=round(mean(storeN[,3]-storeN[,2],na.rm=TRUE),2)
  Nfail=sims-den
  Nstats=data.frame(MeanNhat=round(Nmean,2),Bias=Nbias,Coverage=Ncover,CIwidth=Nwidth,
                    simFails=Nfail)

  #Calculate W stats
   CapStats=rbind(colMeans(storeCapsTrue),colMeans(storeCapsObs),colMeans(storeCapsTrue)-colMeans(storeCapsObs))
   rownames(CapStats)=c("True Capture Events","Observed Capture Events","NMissing")
   CapStats=cbind(CapStats,rowSums(CapStats))
   colnames(CapStats)=c(paste("T",1:occ,sep=""),"Total")
   CapStats=rbind(CapStats,CapStats[2,]/CapStats[1,])
   rownames(CapStats)[4]="%remain"

  #Calculate hair stats
  HairStats=rbind(colMeans(storeS),colMeans(storeU),colMeans(storeR))
  HairStats=cbind(HairStats,rowSums(HairStats))
  colnames(HairStats)=c(paste("T",1:occ,sep=""),"Total")
  HairStats=rbind(HairStats,HairStats[2,]/HairStats[1,],HairStats[3,]/HairStats[1,])
  rownames(HairStats)=c("Sdotj","Udotj","Rdotj","%remainU","%remainR")
  HairStats=round(HairStats,2)


  return(list(Nstats=Nstats,CapStats=CapStats,HairStats=HairStats,storeN=storeN,
              storeS=storeS,storeU=storeU))
}
