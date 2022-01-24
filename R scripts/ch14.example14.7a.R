## compare marginal models versus conditional models


createData<-function(nc=NC,m=M,beta=BETA,sigma=SIGMA,mu=MU){
  cluster.effects<- rnorm(nc,sd=sigma)
  gamma<- rep(cluster.effects,each=m)
  if (m/2 - floor(m/2) > .1) stop("m must be an even positive integer")
  trt<-rep(c(rep(0,m/2),rep(1,m/2)),nc)
  linpred<- mu + gamma + beta*trt
  y<-rbinom(nc*m,1,plogis(linpred))
  list(gamma=gamma, trt=trt, y=y, cid=rep(1:nc,each=m))
}

#d<-createData(nc=500,m=2,beta=1,sigma=1,mu=-1)
#d<-as.data.frame(d)
#d$cid<- factor(d$cid)
#gout<-glm(y~cid+trt,family=binomial, data=d)
#t0<-proc.time()
#gout<-glmfunc(d)
#t1<-proc.time()
#system.time(gout2<-glmfunc2(d) )
#t2<-proc.time()
#t1-t0
#t2-t1
#smGLM(NSIM=100,NC=500,M=2,BETA=1,MU=-1,SIGMA=1)

#d<-as.data.frame(d)

library(speedglm)
glmfunc2<-function(d){
  d<-as.data.frame(d)
  gout<-glm(y~factor(cid)+trt,family=binomial, data=d)
  #gci<-confint(gout)["trt",]
  stderr<-sqrt(vcov(gout)["trt","trt"])

  gci<-coef(gout)["trt"] + c(-1,1)*qnorm(.975)*stderr
  
  c(est=coef(gout)["trt"],lo=gci[1],hi=gci[2],p.value=summary(gout)$coef["trt","Pr(>|z|)"])
}





library(lme4)
glmmixfunc<-function(d){
  greout<-glmer(y~1+trt+(1|cid),data=d,family=binomial())
  sdbeta<- sqrt(vcov(greout)["trt","trt"])
  beta<- fixef(greout)["trt"]
  z<- qnorm(.975)
  lo<- beta - z*sdbeta
  hi<- beta + z*sdbeta
  p.value<- 2*(1-pnorm(abs(beta/sdbeta)))
  #gci<-confint(greout)["trt",]
  c(est=beta,lo=lo,hi=hi,p.value=p.value)
}

library(saws)
library(survival)
clogitfunc<-function(d){
  D<-as.data.frame(d)
  cout<-clogit(y~trt+strata(cid), data=D)
  ci<-confint(cout)
  c(est=coef(cout),lo=ci[1],hi=ci[2],p.value=summary(greout)$coef["trt","Pr(>|z|)"])
}

#clogitfunc(d)


library(gee)
library(saws)
sawsfunc<-function(d,...){
  geeout<- mgee(y~trt,id=cid,family=binomial,data=as.data.frame(d),...)
  sout<-saws(geeout)
  c(est=sout$coef[2,1],lo=sout$conf.int[2,1],hi=sout$conf.int[2,2],p.value=sout$p.value[2])
}  


## Appear to get the same answers with independence and exchangeable.
#d<-createData(m=6)
#sawsfunc(d)
#sawsfunc(d,corstr="exchangeable")

#gee(y~trt,id=cid,family=binomial,data=as.data.frame(d),corstr = "exchangeable")
#gee(y~trt,id=cid,family=binomial,data=as.data.frame(d),corstr = "independence")




simGLM<-function(NSIM=10,NC=100,M=2,BETA=1,MU=-1,SIGMA=1){
  out<- matrix(NA,NSIM,4)
  for (i in 1:NSIM){
    d<-createData(nc=NC,m=M,beta=BETA,sigma=SIGMA,mu=MU)
    out[i,]<-glmfunc(d) 
  }
  errhi<-sum( out[,3]<BETA)/NSIM
  errlo<- sum( out[,2]>BETA)/NSIM
  coverage<- sum( out[,2]<=BETA & out[,3]>=BETA )/NSIM
  mean<- mean(out[,1])
  power<- sum(out[,4]<=0.05)/NSIM
  qvals<- quantile(out[,1],probs=c(.25,.5,.75))
  attributes(qvals)<-NULL
  list(out=out, summary=c(
    errhi=errhi,errlo=errlo,coverage=coverage,
    mean=mean,power=power,q25=qvals[1],q50=qvals[2],q75=qvals[3])
  )
}

# Example 14.7a
nsim<-10^4
set.seed(10201)
s2<-simGLM(NSIM=nsim,NC=100,M=2,BETA=1,MU=-1,SIGMA=1)
s2$summary
set.seed(4813)
s4<-simGLM(NSIM=nsim,NC=100,M=4,BETA=1,MU=-1,SIGMA=1)
s4$summary
set.seed(685911)
s6<-simGLM(NSIM=nsim,NC=100,M=6,BETA=1,MU=-1,SIGMA=1)
s6$summary
set.seed(594921)
s10<-simGLM(NSIM=nsim,NC=100,M=10,BETA=1,MU=-1,SIGMA=1)
s10$summary
set.seed(499911)
s20<-simGLM(NSIM=nsim,NC=100,M=20,BETA=1,MU=-1,SIGMA=1)
s20$summary

t0<-proc.time()
# to do 1000, it should take about 7 hours
nsim<-1000
set.seed(104011)
system.time(s2.1000<-simGLM(NSIM=nsim,NC=1000,M=2,BETA=1,MU=-1,SIGMA=1) )
s2.1000$summary
proc.time() - t0


#save(s2,s4,s6,s10,s20,s2.1000,file="chModMarginalvsConditionalResults1.Rdata")


###########################################################################
#  Second Simulation
###########################################################################

set.seed(040211)
NSIM<-10000
r.glm<- r.glmmix<-r.saws<-r.clogit<- matrix(NA,NSIM,4)
for (i in 1:NSIM){
  d<-createData(nc=100,m=2,beta=1,sigma=1,mu=-1)
  r.glm[i,]<- glmfunc(d) 
  r.glmmix[i,]<-glmmixfunc(d)
  r.clogit[i,]<-clogitfunc(d)
  r.saws[i,]<- sawsfunc(d)
  print(paste("i=",i))
}

## 

save(r.glm,r.glmmix,r.clogit,r.saws,file="chModMarginalvsConditionalResults2.Rdata")

quantile(r.glm[,1],probs=c(1:9/10))
quantile(r.glm[,1],probs=c(.25,.75))
quantile(r.clogit[,1],probs=c(.25,.5,.75))
quantile(r.saws[,1],probs=c(1:9/10))
quantile(r.glmmix[,1],probs=c(.25,.5,.75))

par(mfrow=c(1,3))
hist(r.glm[,1],xlim=c(-2,4))
hist(r.clogit[,1],xlim=c(-2,4))
hist(r.saws[,1],xlim=c(-2,4))

mean(r.glm[,1])
mean(r.glmmix[,1])
mean(r.clogit[,1])
mean(r.saws[,1])

k<- (16*sqrt(3)/(15*pi))
k^2
1/sqrt(1+ k^2)
1/sqrt(1+0.346)



### did not use tabSim...need to edit it if want to use it for a table

tabSim<-function(sresult,aslog=TRUE){
  
  nsim<-nrow(sresult$ratio)
  if (aslog){
    ratio<-log10(sresult$ratio)
    lo<-log10(sresult$lo)
    hi<-log10(sresult$hi)
    RATIO<- log10(sresult$trueRatio)
    
  } else {
    ratio<-sresult$ratio
    lo<-sresult$lo
    hi<-sresult$hi
    RATIO<-sresult$trueRatio
    
  }
  
  
  out<-matrix(NA,2,9,dimnames=list(dimnames(ratio)[[2]],c("Pct Convg","bias","var","mse","cover","mean CI len","med CI len",
                                                          "nMissEst","nMissCI")))
  
  
  isMiss<- function(x){ is.na(x) | x==Inf | x==-Inf }
  nmiss<-function(x){ length(x[isMiss(x)])}
  
  for (j in 1:2){
    rj<- ratio[,j]
    out[j,1]<- 100 - 100*(nmiss(rj)/nsim)
    I<-  !isMiss(rj)
    out[j,2]<- mean(rj[I] - RATIO )
    out[j,3]<- ((nsim-1)/nsim)*var(rj[I])
    out[j,4]<-  mean( (rj[I]-RATIO)^2 )
    # nnm=n not missing
    nnm<- length(rj[I])
    out[j,5]<-  100* ( length(rj[I][lo[I,j]<=RATIO & hi[I,j]>=RATIO])/nnm    ) 
    J<- !isMiss(lo[,j]) & !isMiss(hi[,j])
    out[j,6]<- mean(hi[J,j] - lo[J,j])
    out[j,7]<- median(hi[J,j] - lo[J,j])
    out[j,8]<- nsim - nnm
    out[j,9]<- nsim- length(rj[J])
  }
  
  
  out
}


