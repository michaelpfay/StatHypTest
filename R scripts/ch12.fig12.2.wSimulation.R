# Figure 12.2

# Redo Simulation of Michael et al 2019 Biometrics DOI: 10.1111/biom.12998
#source("metaAnalysis.R")
library(asht)
#library(meta)
library(rma.exact)

set.seed(10301)
NSIM<- 10^4
Kvals<- 3:20

t0<-proc.time()

## define functions need for simulation
##
createData<-function(nsim,K, tau0sq=0, mu0=0){
  sigk<- 1 + 4*(1:K -1)/(K-1)
  out<-matrix(NA,nsim,K)
  SD<- sqrt(tau0sq+sigk^2)
  for (i in 1:nsim){
    out[i,]<- rnorm(K,mean=mu0, sd=SD )
  }
  list(y=out,sigk=sigk,tau0sq=tau0sq)
}

dorow<-function(y,v){
  k<-length(y)
  CI<-matrix(NA,6,2,dimnames=list(mnames,c("lo","hi")))
  CI[1,]<-metaNorm(y,v,method="fixed")$conf.int
  CI[2,]<-metaNorm(y,v,method="DL",df=Inf)$conf.int
  CI[3,]<-metaNorm(y,v,method="DL",df=k-1)$conf.int
  CI[4,]<-metaNorm(y,v,method="PM",df=Inf)$conf.int
  CI[5,]<-metaNorm(y,v,method="PM",df=k-1)$conf.int
  #CI[6,]<-rma.exact.fast(y,v, plot=FALSE)[1,]
  list(lo=CI[,1],hi=CI[,2]) 
}


propReject<-function(x,b0=0,direction="less"){ 
  if (direction=="greater"){
    out<-length(x[x>=b0])/length(x)
  } else {
    out<-length(x[x<=b0])/length(x)
  }
  out
}



## SIMULATION 
doSimulation<-function(tau2){
  mnames<-c("fixed","DL norm","DL t","PM norm","PM t","rma exact")
  
  errlo<-errhi<- cilen<-matrix(NA,length(Kvals),length(mnames),dimnames=list(Kvals,mnames))
  
  for (j in 1:length(Kvals)){
    print(paste("K=",Kvals[j], " tau2=",tau2))
    x<-createData(NSIM,Kvals[j], tau0sq = tau2)
    #rma.exact.fast(x$y[1,],x$sigk^2)
    
    lower<-upper<-matrix(NA,NSIM,6,dimnames=list(NULL,mnames))
    for (i in 1:NSIM){
      temp<-dorow(x$y[i,],x$sigk^2)
      lower[i,]<-temp$lo
      upper[i,]<-temp$hi
    }
    cilen[j,]<- apply(upper - lower,2,mean)
    errlo[j,]<-apply(lower,2,propReject,direction="greater")
    errhi[j,]<-apply(upper,2,propReject,direction="less")
  }
  list(cilen=cilen,errlo=errlo,errhi=errhi)
}


simResult.tau2_0<- doSimulation(0)
simResult.tau2_0
#save(simResult.tau2_0,file="chStr.simMetaAnalysis.Tau2.0.Rdata")
simResult.tau2_12.5<- doSimulation(12.5)
simResult.tau2_12.5
#save(simResult.tau2_12.5,file="chStr.simMetaAnalysis.Tau2.12.5.Rdata")
simResult.tau2_25<- doSimulation(25)
simResult.tau2_25
#save(simResult.tau2_25,file="chStr.simMetaAnalysis.Tau2.25.Rdata")


#load("chStr.simMetaAnalysis.Tau2.0.Rdata")
#load("chStr.simMetaAnalysis.Tau2.12.5.Rdata")
#load("chStr.simMetaAnalysis.Tau2.25.Rdata")





plotSim<-function(sim,TITLE="",dofix=TRUE){
  err<- sim$errlo + sim$errhi
  coverage<- 1 - err
  plot(c(3,20),c(.75,1),type="n",xlab="K",ylab="Coverage", main=TITLE)
  lines(3:20,coverage[,"DL norm"],lwd=3,lty=1,col="gray")
  lines(3:20,coverage[,"DL t"],lwd=6,lty=1,col="gray")
  lines(3:20,coverage[,"PM t"],lwd=3,lty=1)
  lines(c(0,30),c(.95,.95),lwd=1,lty=3)
  if (dofix){
    lines(3:20,coverage[,"fixed"],lwd=3,lty=3,col="gray")
  }
}

par(mfrow=c(1,3))
plotSim(simResult.tau2_0,expression(tau^2 == 0))
legend("bottomleft", legend=c("PM t","DL t","DL norm","fixed"),lwd=c(3,6,3,3),lty=c(1,1,1,3),col=c("black","gray","gray","gray"))
plotSim(simResult.tau2_12.5,expression(tau^2 == 12.5))
plotSim(simResult.tau2_25,expression(tau^2 == 25))

#dev.print(pdf,file="../book/graph/chStrSimMetaAnalysis.pdf")

proc.time() - t0