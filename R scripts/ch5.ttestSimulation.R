# Run simulation for Figures 5.2 and 5.3
# Might need some tweeking to reproduce, since I ran it in two batches
# and may have saved results in files with different names

#setwd("S:/BRB/Staff Folders/Mike Fay/book/ASHT/R")


## Feb/13/2015: Was going to change log gamma to have 
## mean at zero, decided not to change, but reran simulation
## for all except bootstrap, since I accidentally erased the
## old simulation before we decided not to change it.
t0<-proc.time()
NSIM<-10^6
library(boot)



ttestSimulation<-function(n,nsim,nulldist,nullmean,altmean,conf.level=0.95,boot=FALSE,btype="perc"){
    ci<-matrix(NA,nsim,2)
    if (boot){
      Mean<-function(x,i){ mean(x[i]) }
      ## for btype="stud" use mean.fun
      mean.fun <- function(x, i){    
        m <- mean(x[i])
        n <- length(i)
        v <- (n-1)*var(x[i])/n^2
        c(m, v)
      }
    } 
    for (i in 1:nsim){
        if (boot){
          if (btype=="perc"){
            ci[i,]<-boot.ci(boot(nulldist(n),Mean,R=1999,stype="i"),type=btype)$percent[4:5]
          } else if (btype=="bca"){
            ci[i,]<-boot.ci(boot(nulldist(n),Mean,R=1999,stype="i"),type=btype)$bca[4:5]
          } else if (btype=="stud"){
            ## Aug/09/2017--added student-t CI because of Owen, 1988 
            ## see articles/empirical likelihood/Owen.1988.pdf
            ci[i,]<-boot.ci(boot(nulldist(n),mean.fun,R=1999,stype="i"),type=btype)$student[4:5]
          }
        } else {
            ci[i,]<-t.test(nulldist(n))$conf.int
        }
    }
    x<-(1:nsim)
    means<-c(nullmean,altmean)
    nm<-length(means)
    rejectLow<-rejectHi<-rep(NA,nm)
    for (i in 1:length(means)){
        rejectLow[i]<- 100*length(x[ci[,1]>means[i]])/nsim
        rejectHi[i]<- 100*length(x[ci[,2]<means[i]])/nsim
    }
    names(rejectLow)<-names(rejectHi)<-round(means,4)
    list(rejectLow=rejectLow,rejectHi=rejectHi)
}

#npois<-function(n,u=3){
#    rpois(n,u)
#}

#x<-ttestSimulation(50,10^4,npois,3,c(2.5,3.5,4,5))
#x


set.seed(103021)
rnorm1<-function(n){ rnorm(n,sd=1/qnorm(.75)) }
sim.norm.20<-ttestSimulation(20,NSIM,rnorm1,0,c(-1,-.5,.5,1))
sim.norm.50<-ttestSimulation(50,NSIM,rnorm1,0,c(-1,-.5,.5,1))
sim.norm.200<-ttestSimulation(200,NSIM,rnorm1,0,c(-1,-.5,.5,1))

# Cauchy
rtdf1<-function(n,df=1){
    rt(n,df)
}
sim.tdf1.20<-ttestSimulation(20,NSIM,rtdf1,0,c(-1,-.5,.5,1))
sim.tdf1.50<-ttestSimulation(50,NSIM,rtdf1,0,c(-1,-.5,.5,1))
sim.tdf1.200<-ttestSimulation(200,NSIM,rtdf1,0,c(-1,-.5,.5,1))

# t df=2
rtdf2<-function(n,df=2){
    rt(n,df)
}
sim.tdf2.20<-ttestSimulation(20,NSIM,rtdf2,0,c(-1,-.5,.5,1))
sim.tdf2.50<-ttestSimulation(50,NSIM,rtdf2,0,c(-1,-.5,.5,1))
sim.tdf2.200<-ttestSimulation(200,NSIM,rtdf2,0,c(-1,-.5,.5,1))

# Log gamma
source("ch5.ttestSimGetLogGammaParms.R")
rLogGamma<-function(n){
    rlgamma(n,location=loc,shape=SHAPE,scale=s)
}
rLogGamma(10)

sim.logGam.20<-ttestSimulation(20,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1))
sim.logGam.50<-ttestSimulation(50,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1))
sim.logGam.200<-ttestSimulation(200,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1))

sim.norm.20
sim.norm.50
sim.norm.200
sim.tdf1.20
sim.tdf1.50
sim.tdf1.200
sim.tdf2.20
sim.tdf2.50
sim.tdf2.200
sim.logGam.20
sim.logGam.50
sim.logGam.200

save.image("simulation1.RData")

NSIM<-10^6
set.seed(103101)
sim.logGam.20.bootperc<-ttestSimulation(20,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1),boot=TRUE)
sim.logGam.20.bootbca<-ttestSimulation(20,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1),boot=TRUE,btype="bca")
set.seed(949381)
t0<-proc.time()
NSIM<-10^5
sim.logGam.20.bootstud<-ttestSimulation(20,NSIM,rLogGamma,logGammaMean,c(-1,-.5,.5,1),boot=TRUE,btype="stud")
t1<-proc.time()
t1-t0

sim.logGam.20.bootperc
sim.logGam.20.bootbca
sim.logGam.20.bootstud
save.image("simulation1boot.RData")
#load("simulation1boot.RData")

proc.time()-t0