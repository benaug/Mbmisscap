 #'This function asseses some frequentist properties when naively fitting Mb to capture histories that are missing data.
#'@param N an integer that is the true population size
#'@param p a vector of size 2 containing the capture and recapture probabilities
#'@param occ an integer that is the number of capture occasions to simulate
#'@param K an integer that is the number of traps to simulate
#'@param traptype a character string specifying the trap type.  "single" 
#'allows individuals to be caught in at most 1 trap per occasion and "multi" allows
#'individuals to be caught in multiple traps per occasion.
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
#'@param sims the number of simulated data sets to generate and fit Mb to
#'@param cores the number of cores to do the simulation on.  A typical for loop is used if cores=1 and
#'the foreach package is used if cores>1
#'@return a list with the simulation results (expand later)
#'@include simMBmisscap
#'@export

assess=function(N,p,occ,K,traptype,btype="global",lambda_h,lambda_c=NULL,delta=NULL,kappa_h=NULL,
                kappa_t=NULL,cluster=FALSE,kappa_c=NULL,alpha=NULL,sims=100,cores=1){
  if(cores==1){
    #preallocate storage structures
    storeN=matrix(NA,nrow=sims,ncol=3)
    storep=rep(NA,sims)
    storec=rep(NA,sims)
    storeFail=rep(NA,sims)
    storeCapsTrue=matrix(NA,nrow=sims,ncol=occ)
    storeCapsObs=matrix(NA,nrow=sims,ncol=occ)
    storeFirstcapstrue=matrix(NA,nrow=sims,ncol=occ)
    storeMiss1caps=matrix(NA,nrow=sims,ncol=occ)
    storeMdottrue=rep(NA,sims)
    storeMdotobs=rep(NA,sims)
    storeMt1true=rep(NA,sims)
    storeMt1obs=rep(NA,sims)
    storeS=matrix(NA,nrow=sims,ncol=occ)
    storeU=matrix(NA,nrow=sims,ncol=occ)
    storeR=matrix(NA,nrow=sims,ncol=occ)
    storepstats=array(NA,dim=c(5,occ,sims))
    for(i in 1:sims){
      data=simMBmisscap(N,p,occ,K,traptype,btype,lambda_h,lambda_c,delta,kappa_h,
                        kappa_t,cluster,kappa_c,alpha)
      Mdata=as.data.frame(apply(data$Wobs2D,1,paste,collapse=""))
      colnames(Mdata)="ch"
      Mdata$ch=as.character(Mdata$ch)
      out=RMark::mark(Mdata,model="Closed",output=FALSE,delete=TRUE) #Mb
      storeN[i,]=unlist(out$results$derived[c(1,3,4)])
      storep[i]=RMark::get.real(out,"p")[1]
      storec[i]=RMark::get.real(out,"c")[1]
      storeFail[i]=data$Wstats$fail
      storeCapsTrue[i,]=colSums(data$W2D)
      storeCapsObs[i,]=colSums(data$Wobs2D)
      storeFirstcapstrue[i,]=data$Wstats$firstcapsTrue
      storeMiss1caps[i,]=data$Wstats$miss1caps
      storeMdottrue[i]=data$Wstats$Mdot_true
      storeMdotobs[i]= data$Wstats$Mdot_obs
      storeMt1true[i]=data$Wstats$Mt1_true
      storeMt1obs[i]=data$Wstats$Mt1_obs
      storeS[i,]=data$HairStats$Sdotj
      storeU[i,]=data$HairStats$Udotj
      storeR[i,]=data$HairStats$Rdotj
      storepstats[,,i]=data$pstats
    }
    den=sims
  }else{
    #load multicore packages
    require(foreach)
    require(snow)
    require(doSNOW)
    require(RMark)
    cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp)
    #clusterExport(cl.tmp, varlist=c("alpha","btype","cluster","delta"), envir = environment())
    clusterExport(cl.tmp, list=ls(environment()), envir = environment())
    
    out2=foreach(i=1:sims,.packages=c("VGAM","RMark"),.export="simMBmisscap",.errorhandling = "remove") %dopar% {
      data=simMBmisscap(N,p,occ,K,traptype,btype,lambda_h,lambda_c,delta,kappa_h,
                        kappa_t,cluster,kappa_c,alpha)
      Mdata=as.data.frame(apply(data$Wobs2D,1,paste,collapse=""))
      colnames(Mdata)="ch"
      Mdata$ch=as.character(Mdata$ch)
      out=mark(Mdata,model="Closed",output=FALSE,delete=FALSE,threads=1) #Mb 
      f0=get.real(out,"f0",se=TRUE)[3:4]
      c=exp(1.96*sqrt((log(1+f0[2]^2/(f0[1]^2)))))
      Nhat=unlist(c(data$n+f0[1],data$n+f0[1]/c,data$n+f0[1]*c))
      #Nhat=out$results$derived[c(1,3,4)]
      params=list(Nhat,
                  get.real(out,"p")[1],
                  get.real(out,"c")[1])
      
      return(list(params,
                  data$Wstats$fail,
                  colSums(data$W2D),
                  colSums(data$Wobs2D),
                  data$Wstats$firstcapsTrue,
                  data$Wstats$miss1caps,
                  data$Wstats$Mdot_true,
                  data$Wstats$Mdot_obs,
                  data$Wstats$Mt1_true,
                  data$Wstats$Mt1_obs,
                  data$HairStats$Sdotj,
                  data$HairStats$Udotj,
                  data$HairStats$Rdotj,
                  data$pstats))
    }
    stopCluster(cl.tmp)
    cleanup(ask=FALSE,lx=ls(.GlobalEnv))
    den=length(out2) #successful fits RMark didn't loose track of
    #preallocate storage structures
    storeN=matrix(NA,nrow=den,ncol=3)
    storep=rep(NA,den)
    storec=rep(NA,den)
    storeFail=rep(NA,den)
    storeCapsTrue=matrix(NA,nrow=den,ncol=occ)
    storeCapsObs=matrix(NA,nrow=den,ncol=occ)
    storeFirstcapstrue=matrix(NA,nrow=den,ncol=occ)
    storeMiss1caps=matrix(NA,nrow=den,ncol=occ)
    storeMdottrue=rep(NA,den)
    storeMdotobs=rep(NA,den)
    storeMt1true=rep(NA,den)
    storeMt1obs=rep(NA,den)
    storeS=matrix(NA,nrow=den,ncol=occ)
    storeU=matrix(NA,nrow=den,ncol=occ)
    storeR=matrix(NA,nrow=den,ncol=occ)
    storepstats=array(NA,dim=c(5,occ,den))
    for(i in 1:den){
      storeN[i,]=out2[[i]][[1]][[1]]
      storep[i]=out2[[i]][[1]][[2]]
      storec[i]=out2[[i]][[1]][[3]]
      storeFail[i]=out2[[i]][[2]]
      storeCapsTrue[i,]=out2[[i]][[3]]
      storeCapsObs[i,]=out2[[i]][[4]]
      storeFirstcapstrue[i,]=out2[[i]][[5]]
      storeMiss1caps[i,]=out2[[i]][[6]]
      storeMdottrue[i]=out2[[i]][[7]]
      storeMdotobs[i]=out2[[i]][[8]]
      storeMt1true[i]=out2[[i]][[9]]
      storeMt1obs[i]=out2[[i]][[10]]
      storeS[i,]=out2[[i]][[11]]
      storeU[i,]=out2[[i]][[12]]
      storeR[i,]=out2[[i]][[13]]
      storepstats[,,i]=out2[[i]][[14]]
    }
  }
  #Calcuate N stats
  Nmean=mean(storeN[,1],na.rm=TRUE)
  dir=Nmean[1]-N
  Nbias=ifelse(dir>0,paste("+",round(100*(Nmean[1]/N-1),2),"%"),paste("-",round(100*(1-Nmean[1]/N),2),"%"))
  Ncover=round(sum(storeN[,2]<N&storeN[,3]>N,na.rm=TRUE)/den,2)
  Nwidth=round(mean(storeN[,3]-storeN[,2],na.rm=TRUE),2)
  Nfail=sims-den
  LikFails0=sum(storeFail<=0)
  LikFails100=sum(storeFail<=100)
  Nstats=data.frame(MeanNhat=round(Nmean,2),Bias=Nbias,Coverage=Ncover,CIwidth=Nwidth,
                    simFails=Nfail,LikFails0=LikFails0,LikFails100=LikFails100)
  
  #Calculate capstats stats
  CapStats=rbind(colMeans(storeCapsTrue),colMeans(storeCapsObs),colMeans(storeCapsTrue)-colMeans(storeCapsObs))
  rownames(CapStats)=c("True Capture Events","Observed Capture Events","NMissing")
  CapStats=cbind(CapStats,rowSums(CapStats))
  colnames(CapStats)=c(paste("T",1:occ,sep=""),"Total")
  CapStats=rbind(CapStats,CapStats[2,]/CapStats[1,])
  rownames(CapStats)[4]="%remain"
  CapStats=rbind(CapStats,c(colMeans(storeFirstcapstrue),sum(colMeans(storeFirstcapstrue))))
  CapStats=rbind(CapStats,c(colMeans(storeFirstcapstrue-storeMiss1caps),sum(colMeans(storeFirstcapstrue-storeMiss1caps))))
  CapStats=rbind(CapStats,CapStats[5,]-CapStats[6,])
  CapStats=rbind(CapStats,CapStats[6,]/CapStats[5,])
  CapStats=round(CapStats,2)
  rownames(CapStats)[5:8]=c("True First Captures","Observed First Captures","NMissing","%remain")
  
  #Calculate hair stats
  HairStats=rbind(colMeans(storeS),colMeans(storeU),colMeans(storeR))
  HairStats=cbind(HairStats,rowSums(HairStats))
  colnames(HairStats)=c(paste("T",1:occ,sep=""),"Total")
  HairStats=rbind(HairStats,HairStats[2,]/HairStats[1,],HairStats[3,]/HairStats[1,])
  rownames(HairStats)=c("Sdotj","Udotj","Rdotj","%remainU","%remainR")
  HairStats=round(HairStats,2)
  
  #calculate p stats
  pstats=apply(storepstats,1:2,mean)
  rownames(pstats)=c("p_true","p_obs","p_obs_true","p_obs_false","percent_false")
  colnames(pstats)=paste("T",1:occ,sep="")
  pstats=round(pstats,2)
  
  #Calculate Mb stats
  MbStats=matrix(c(mean(storeMt1true),mean(storeMt1obs),mean(storeMdottrue),mean(storeMdotobs)),nrow=2,byrow=TRUE)
  rownames(MbStats)=c("M_t+1","M.")
  colnames(MbStats)=c("True","Observed")
  return(list(Nstats=Nstats,CapStats=CapStats,HairStats=HairStats,MbStats=MbStats,pstats=pstats,storeN=storeN,storep=storep,storec=storec,storeFail=storeFail,
              storeMdottrue=storeMdottrue,storeMdotobs=storeMdotobs,storeMt1true=storeMt1true,
              storeMt1obs=storeMt1obs,storeS=storeS,storeU=storeU))
}