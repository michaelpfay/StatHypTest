# examples 14.10, 14.11, 14.12
# Simulate model building using stepwise
n<-100
nsim<-1E4
t0<-proc.time()
set.seed(201019)
library(mvtnorm)

k<-11
rho<- 0

getpA<-function(d){
  k<-ncol(d)
  D<- data.frame(d)
  pval<-rep(NA,k-1)
  for (i in 2:k){
    f<-formula(paste0("X1~X",i))
    lout<-lm(f,data=D)
    pval[i-1]<-summary(lout)$coef[2,"Pr(>|t|)"]
  }
  p.value<-min(pval)
  p.value
}

getpB<-function(d){
  k<-ncol(d)
  D<-data.frame(d)
  f<-formula(paste0("X1~",paste0("X",2:k,collapse = "+")))
  lout<-lm(f,data=D)
  pval<-summary(lout)$coef[-1,"Pr(>|t|)"]
  p.value<-min(pval)
  p.value
}
getpC<-function(d){
  k<-ncol(d)
  D<-data.frame(d)
  #l0<-lm(X1~X2,data=D)
  l0<-lm(X1~1,data=D)
  f<-formula(paste0("X1~",paste0("X",2:k,collapse = "+")))
  s<-step(l0,scope=list(lower=X1~1,upper=f),trace=0)
  pval<-summary(s)$coef[-1,"Pr(>|t|)"]
  if (is.null(pval)){
    p.value<-1
  }  else {
    p.value<-min(pval)
  }
  p.value
}

pa<-pb<-pc<-rep(NA,nsim)

for (i in 1:nsim){
  d<-rmvnorm(n,sigma=(1-rho)*diag(k) + matrix(rho,k,k) )
  pa[i]<-getpA(d)
  pb[i]<-getpB(d)
  pc[i]<-getpC(d)
   print(paste("i=",i))
}

power<-function(p,alpha=0.05){ length(p[p<=alpha])/length(p)}

power(pa)
power(pb)
power(pc)

proc.time()-t0


##########################################################
# Simulation 2: 1 covariate of interest, the other 9 
#    are nuisance covariates
##########################################################

set.seed(409201)
n<-100
nsim<-1E4
t0<-proc.time()
library(mvtnorm)

k<-11
rho<- 0



getpA2<-function(d){
  x1<-d[,1]
  x2<-d[,2]
  lout<-lm(x1~x2)
  p.value<- summary(lout)$coef[2,"Pr(>|t|)"]
  p.value
}

getpB2<-function(d){
  x1<-d[,1]
  x2<-d[,2]
  l0<-lm(x1~x2)
  pval<-rep(NA,ncol(d)-1)
  pval[1]<- summary(l0)$coef[2,"Pr(>|t|)"]
  for (j in 3:ncol(d)){
    xj<- d[,j]
    l0<-lm(x1~x2+xj)
    pval[j-1]<- summary(l0)$coef[2,"Pr(>|t|)"]
  }
  p.value<-min(pval)
  p.value
}
getpC2<-function(d){
  k<-ncol(d)
  D<-data.frame(d)
  l0<-lm(X1~X2,data=D)
  f<-formula(paste0("X1~",paste0("X",2:k,collapse = "+")))
  s<-step(l0,scope=list(lower=X1~X2,upper=f),trace=0)
  p.value<-summary(s)$coef["X2","Pr(>|t|)"]
  p.value
}




pa2<-pb2<-pc2<-rep(NA,nsim)

for (i in 1:nsim){
  d<-rmvnorm(n,sigma=(1-rho)*diag(k) + matrix(rho,k,k) )
  pa2[i]<-getpA2(d)
  pb2[i]<-getpB2(d)
  pc2[i]<-getpC2(d)
}

power<-function(p,alpha=0.05){ length(p[p<=alpha])/length(p)}

power(pa2)
power(pb2)
power(pc2)

proc.time()-t0




##########################################################
# Simulation 3: 1 covariate of interest, the other 9 
#    are nuisance covariates, but correlation among nuisance
##########################################################

set.seed(583921)
n<-100
nsim<-1E4
t0<-proc.time()
library(mvtnorm)

k<-11
rho<- 0.5


pa3<-pb3<-pc3<-rep(NA,nsim)

for (i in 1:nsim){
  Sigma<- (1-rho)*diag(k) + matrix(rho,k,k)
  Sigma[2,1]<-Sigma[1,2]<-0
  d<-rmvnorm(n,sigma=Sigma)
  pa3[i]<-getpA2(d)
  pb3[i]<-getpB2(d)
  pc3[i]<-getpC2(d)
}

power<-function(p,alpha=0.05){ length(p[p<=alpha])/length(p)}

power(pa3)
power(pb3)
power(pc3)

proc.time()-t0

#save(pa,pb,pc,pa2,pb2,pc2,pa3,pb3,pc3,file="chModStepwiseSimResults.RData")

