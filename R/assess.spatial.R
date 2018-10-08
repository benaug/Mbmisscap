#'This function asseses some frequentist properties when naively fitting Mb to spatial capture histories that are missing captures due to subsampling and/or failed DNA amplification.
#'Models are fit via the secr package.
#'#'@param N an integer that is the true population size
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
#'@export

assess.spatial=function(N,occ,lambda,sigma,locs,buff,spacing,btype="global",lambda_h,lambda_c=NULL,
                        delta=NULL,kappa_h=NULL,kappa_t=NULL,cluster=FALSE,kappa_c=NULL,
                        alpha=NULL,sims=100,cores=1){
  #Get starting values for betas
  #Choose model
  if(btype=="trap"){
    mod=g0~bk
    D=(N/prod(apply(locs,2,function(x){max(x)-min(x)})))*1000
    starts=c(log(D),qlogis(lambda[1]),qlogis(lambda[2])-qlogis(lambda[1]),log(sigma))
  }else if(btype=="global"){
    mod=g0~b
    D=(N/prod(apply(locs,2,function(x){max(x)-min(x)})))*1000
    starts=c(log(D),qlogis(lambda[1]),qlogis(lambda[2])-qlogis(lambda[1]),log(sigma))
  }else if(btype=="none"){
    mod=g0~1
    D=(N/prod(apply(locs,2,function(x){max(x)-min(x)})))*1000
    starts=c(log(D),qlogis(lambda[1]),log(sigma))
  }
  if(cores==1){
    #preallocate storage structures
    storeN=matrix(NA,nrow=sims,ncol=3)
    storeCapsTrue=matrix(NA,nrow=sims,ncol=occ)
    storeCapsObs=matrix(NA,nrow=sims,ncol=occ)
    storeFirstcapstrue=matrix(NA,nrow=sims,ncol=occ)
    storeMiss1caps=matrix(NA,nrow=sims,ncol=occ)
    storeS=matrix(NA,nrow=sims,ncol=occ)
    storeU=matrix(NA,nrow=sims,ncol=occ)
    storeR=matrix(NA,nrow=sims,ncol=occ)
   # storepstats=array(NA,dim=c(5,occ,sims))
    for(i in 1:sims){
      data=simSpatial(N,occ,lambda,sigma,locs,buff[1],btype,lambda_h,lambda_c=lambda_c,
                        delta=delta,kappa_h=kappa_h,kappa_t=kappa_t,cluster=cluster,
                      kappa_c=kappa_c,alpha=alpha)
      n=sum(data$Wobs)
      input=data.frame(session=rep("a",n),ID=rep(0,n),Occasion=rep(0,n),Detector=rep(0,n))
      caps=which(data$Wobs>0,arr.ind=TRUE)
      for(j in 1:nrow(caps)){
        input[j,c(2,3,4)]=c(caps[j,1],caps[j,2],caps[j,3])
      }
      traps=data.frame(num=1:nrow(locs),x=locs[,1],y=locs[,2])
      traps=read.traps(data=traps,detector="proximity")
      capthist=make.capthist(input,traps)
      #plot(capthist,tracks=TRUE)
      ## generate habitat mask
      fitmask=make.mask (traps, buffer = buff[2], spacing = spacing)
      secr.model=secr.fit(capthist,model=mod,mask=fitmask,start=starts)
      predmask=make.mask (traps, buffer = buff[1], spacing = spacing)
      Nhat=region.N(secr.model,region=predmask)
      RN=as.numeric((Nhat[2,c(1,3,4)]))
      storeN[i,]=RN
      storeCapsTrue[i,]=colSums(data$W2D)
      storeCapsObs[i,]=colSums(data$Wobs2D)
      storeFirstcapstrue[i,]=data$Wstats$firstcapsTrue
      storeMiss1caps[i,]=data$Wstats$miss1caps
      storeS[i,]=data$HairStats$Sdotj
      storeU[i,]=data$HairStats$Udotj
      storeR[i,]=data$HairStats$Rdotj
    }
  }else{
    #load multicore packages
    require(foreach)
    require(snow)
    require(doSNOW)
    require(secr)
    cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp)
    #clusterExport(cl.tmp, varlist=c("alpha","btype","cluster","delta"), envir = environment())
    clusterExport(cl.tmp, list=ls(environment()), envir = environment())

    out2=foreach(i=1:sims,.packages=c("VGAM","secr","Mbmisscap","sp","abind"),.export="simSpatial") %dopar% {
      data=simSpatial(N,occ,lambda,sigma,locs,buff[1],btype,lambda_h,lambda_c=lambda_c,
                      delta=delta,kappa_h=kappa_h,kappa_t=kappa_t,cluster=cluster,
                      kappa_c=kappa_c,alpha=alpha)
      n=sum(data$Wobs)
      input=data.frame(session=rep("a",n),ID=rep(0,n),Occasion=rep(0,n),Detector=rep(0,n))
      caps=which(data$Wobs>0,arr.ind=TRUE)
      for(j in 1:nrow(caps)){
        input[j,c(2,3,4)]=c(caps[j,1],caps[j,2],caps[j,3])
      }
      traps=data.frame(num=1:nrow(locs),x=locs[,1],y=locs[,2])
      traps=read.traps(data=traps,detector="proximity")
      capthist=make.capthist(input,traps)
      #plot(capthist,tracks=TRUE)
      ## generate habitat mask
      fitmask=make.mask (traps, buffer = buff[2], spacing = spacing)
      secr.model=secr.fit(capthist,model=mod,mask=fitmask,start=starts)
      predmask=make.mask (traps, buffer = buff[1], spacing = spacing)
      Nhat=region.N(secr.model,region=predmask)
      RN=as.numeric((Nhat[2,c(1,3,4)]))

      return(list(RN,
                  colSums(data$W2D),
                  colSums(data$Wobs2D),
                  data$Wstats$firstcapsTrue,
                  data$Wstats$miss1caps,
                  data$HairStats$Sdotj,
                  data$HairStats$Udotj,
                  data$HairStats$Rdotj))
    }
    stopCluster(cl.tmp)
    #preallocate storage structures
    storeN=matrix(NA,nrow=sims,ncol=3)
    storeCapsTrue=matrix(NA,nrow=sims,ncol=occ)
    storeCapsObs=matrix(NA,nrow=sims,ncol=occ)
    storeFirstcapstrue=matrix(NA,nrow=sims,ncol=occ)
    storeMiss1caps=matrix(NA,nrow=sims,ncol=occ)
    storeS=matrix(NA,nrow=sims,ncol=occ)
    storeU=matrix(NA,nrow=sims,ncol=occ)
    storeR=matrix(NA,nrow=sims,ncol=occ)
    for(i in 1:sims){
      storeN[i,]=out2[[i]][[1]]
      storeCapsTrue[i,]=out2[[i]][[2]]
      storeCapsObs[i,]=out2[[i]][[3]]
      storeFirstcapstrue[i,]=out2[[i]][[4]]
      storeMiss1caps[i,]=out2[[i]][[5]]
      storeS[i,]=out2[[i]][[6]]
      storeU[i,]=out2[[i]][[7]]
      storeR[i,]=out2[[i]][[8]]
    }
  }
  #Calcuate N stats
  Nmean=mean(storeN[,1],na.rm=TRUE)
  dir=Nmean[1]-N
  Nbias=ifelse(dir>0,paste("+",round(100*(Nmean[1]/N-1),2),"%"),paste("-",round(100*(1-Nmean[1]/N),2),"%"))
  Ncover=round(sum(storeN[,2]<N&storeN[,3]>N,na.rm=TRUE)/sims,2)
  Nwidth=round(mean(storeN[,3]-storeN[,2],na.rm=TRUE),2)
  Nstats=data.frame(MeanNhat=round(Nmean,2),Bias=Nbias,Coverage=Ncover,CIwidth=Nwidth)

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

  return(list(Nstats=Nstats,CapStats=CapStats,HairStats=HairStats,storeN=storeN))
}
