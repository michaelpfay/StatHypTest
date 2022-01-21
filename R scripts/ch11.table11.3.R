# Table 11.3


## k-sample pairwise simulation
NSIM<-10^4

library(multcomp)
t0<-proc.time()
library("multcomp")
#source("asht/R/anovaOneWay.R")
#source("asht/R/tukeyWelsch.R")
# M.Fay Jan 2022: simulation was run on an old version of the functions, 
# I think the function arguments and outputs are the same as in 
# the asht R package, but I did not rerun the simulations 
# to check
library(asht)
createData<-function(n,u,ry=rnorm){
  k<-length(u)
  y<-rep(NA,sum(n))
  g<-rep(1:k,times=n)
  cnt<-1
  for (i in 1:k){
    
    y[cnt:(cnt+n[i]-1)]<-ry(n[i],u[i])
    cnt<-cnt+n[i]
  }
  list(y=y,g=factor(g))
}
#d<-createData(c(10^3,10^3,10^3,10^3),c(0,1,2,3)*.01)

get12pvalues<-function(d){
  # get p-value comparing groups 1 and 2
  pout<-rep(NA,5)
  pout[1]<-pairwise.t.test(d$y,d$g,pool.sd=FALSE,p.adjust.method="none")$p.value[1,1]
  pout[2]<-pairwise.t.test(d$y,d$g,pool.sd=FALSE,p.adjust.method="holm")$p.value[1,1]
  pout[3]<-TukeyHSD(aov(d$y~d$g))[[1]][,"p adj"][1]
  pout[4]<-tukeyWelsch(d$y,d$g,method="sr")$pairwise.pvalues["1 and 2"]
  pout[5]<-tukeyWelsch(d$y,d$g,method="aov")$pairwise.pvalues["1 and 2"]
  pout
}



simInput<-list(
  sim1=list(nsim=NSIM,n=c(10,10,10),u=c(1,0,0),rfunc=rnorm),
  sim2=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(1,0,0,0,0,0),rfunc=rnorm),
  sim3=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(1,0,50,50,50,50),rfunc=rnorm),
  sim4=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(0,0,50,50,50,50),rfunc=rnorm),
  sim5=list(nsim=NSIM,n=c(10,10,10),u=c(1,2,2),rfunc=rpois),
  sim6=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(4,2,1.5,1.5,1.5,1.5),rfunc=rpois),
  sim7=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(4,4,1.5,1.5,1.5,1.5),rfunc=rpois),
  sim8=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(4,2,50,50,50,50),rfunc=rpois),
  sim9=list(nsim=NSIM,n=c(8,12,10,8,12,10),u=c(10,10,10,10,10,10),rfunc=rpois)
  )



set.seed(0234501)

dosim<-function(x,alpha=0.05){
  p<-matrix(NA,x$nsim,5)
  for (i in 1:x$nsim){
    d<-createData(x$n,x$u,x$rfunc)
    p[i,]<-get12pvalues(d)
  }
  getPower<-function(pvals){ length(pvals[pvals<=alpha])/x$nsim }
  apply(p,2,getPower)
}


simPower<-matrix(NA,length(simInput),5)

for (i in 1:length(simInput)){
  simPower[i,]<- dosim(simInput[[i]])
  print(paste("sim ",i))
}
dimnames(simPower)<- list(paste0("sim",1:length(simInput)),c("t-test","Holm adj t-test","TukeyHSD","T-W sr","T-W aov"))

library(xtable)

print(xtable(100*simPower,digits=1))

#save(simPower,simInput,file="chKS.ksample.pairwise.Sim.results.RData")
#load("simPower.RData")
t1<-proc.time()
t1-t0
