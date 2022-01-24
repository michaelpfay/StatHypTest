# Section 18.4
# Simulation on Fisher's exact test

library(gsDesign)


creatTableOF<-function(overallTwoSidedAlpha=0.05,K=k,oneSided=TRUE){
  cOF<-rep(NA,length(K))
  ALPHA<-overallTwoSidedAlpha
  nomAlpha<-matrix(NA,length(K),max(K))
  for (i in 1:length(K)){
    g<-gsDesign(k=K[i],test.type=2, alpha=ALPHA/2,sfu="OF")
    cOF[i]<- g$upper$bound[K[i]]
    t<-  1:K[i]/K[i]
    Zbound<- cOF[i]/sqrt(t)
    if (oneSided) nomAlpha[i,1:K[i]]<- (1-pnorm(Zbound))
    else nomAlpha[i,1:K[i]]<- 2*(1-pnorm(Zbound))
  }
  
  #data.frame(K,cOF,nomAlpha)
  nomAlpha
}

nomAlpha<-creatTableOF(K=k)


rejectFunc<-function(d,Alpha=nomAlpha){
  k<- length(d$x1)
  reject<- FALSE
  p<-rep(NA,k)
  for (j in 1:k){
    X1<- sum(d$x1[1:j])
    N1<- sum(d$n1[1:j])
    X2<- sum(d$x2[1:j])
    N2<- sum(d$n2[1:j])
    xmat<- matrix(c(X1,N1-X1,X2,N2-X2),2,2)
    #print(xmat)
    p[j]<-fisher.test(xmat,conf.int = FALSE, alternative="less")$p.value
    reject<- p[j]<=Alpha[j]
    if (reject) break()
  }
  #list(p=p,reject=reject)  
  reject
}

createData<-function(n1=m,n2=m,p1=P1,p2=P2,k=K){
  x1<- rbinom(k,n1,p1)
  x2<-rbinom(k,n2,p2)
  list(x1=x1,x2=x2,n1=rep(n1,k),n2=rep(n2,k))
}


# check 
#d<-createData(p1=.2,p2=.6)
#rejectFunc(d)

t0<-proc.time()
Nsim<- 1e4
NLoop<- 100
M<- 5
K<- 5
set.seed(4020941)


simPower<-function(nsim=Nsim,m=M,k=K,p=P,alphaVector=nomAlpha){
# to avoid recurrive k calls, rename k  
kvalue<-k
reject<-rep(NA,nsim)
for (i in 1:nsim){
  d<-createData(n1=m,n2=m,p1=P[1],p2=P[2],k=kvalue)
  reject[i]<- rejectFunc(d,Alpha=alphaVector)
}

sum(reject)/nsim
}


PNullVals<-c(1:100/101)
PNullVals<-runif(NLoop)

power<-rep(NA,length(PNullVals))

for (h in 1:length(PNullVals)){
  power[h]<- simPower(nsim=Nsim,m=M,k=K,p=rep(PNullVals[h],2),alphaVector=nomAlpha)
}


proc.time()-t0

### Change m and save power after each run...
#power.m.equal.10<-power
#power.m.equal.5<-power


#save(power.m.equal.10,power.m.equal.5,file="chGS.FETsim.Rdata")

#load("chGS.FETsim.Rdata")


range(power.m.equal.10)
range(power.m.equal.5)

