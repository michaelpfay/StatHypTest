# Section 20.3
# Sample size for two-sample binary
library(exact2x2)

ss2armBin<-function(p0,p1,alpha,power,r=1){
  # calculate several sample sizes for 1 sided test of proportions
  
  ssNorm<-function(delta,tau0,tau1,alpha,gamma){
    (tau0*qnorm(1-alpha)+ tau1*qnorm(1-gamma))^2/delta^2
  }
  # difference in proportions
  delta<- p1-p0
  pbar<- (1/(1+1/r))*(p0 + p1/r)
  tau0<-sqrt( (1+1/r)*pbar*(1-pbar))
  tau1<-sqrt( p0*(1-p0) +p1*(1-p1)/r )
  n.diff<- ssNorm(delta,tau0,tau1,alpha,1-power)
  
  # arcsin transformation
  delta<- 2*asin(sqrt(p1)) - 2*asin(sqrt(p0))
  tau0<- tau1<- sqrt( (1+1/r)/r )
  n.asin<-ssNorm(delta,tau0,tau1,alpha,1-power)
  
  list(n.diff=n.diff,n.asin=n.asin)
  
}

ss36<-ss2armBin(.3,.6,.025,.8)

pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=FALSE)$p.value
}

ss36fet<-SS2x2(.3,.6,.025,pvalFunc,power=.80,n1start=ceiling(ss36$n.diff),increaseby = 1)

# Mid-p version
pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=TRUE)$p.value
}

ss36fetmid<-SS2x2(.3,.6,.025,pvalFunc,power=.80,n1start=ceiling(ss36$n.diff),increaseby = 1)

ss15<-ss2armBin(.01,.15,.025,.8)

pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=FALSE)$p.value
}

ss15fet<-SS2x2(.01,.15,.025,pvalFunc,power=.80,n1start=ceiling(ss15$n.diff),increaseby = 1)

# Mid-p version
pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=TRUE)$p.value
}

ss15fetmid<-SS2x2(.01,.15,.025,pvalFunc,power=.80,n1start=ceiling(ss15$n.asin),increaseby = 1)


ss14<-ss2armBin(.1,.4,.025,.8)


pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=FALSE)$p.value
}

ss14fet<-SS2x2(.1,.4,.025,pvalFunc,power=.80,n1start=ceiling(ss14$n.asin),increaseby = 1)

# Mid-p version
pvalFunc<-function(x1,n1,x2,n2){
  fisher.exact(cbind(c(x1,n1-x1),c(x2,n2-x2)),conf.int=FALSE,alternative="less",midp=TRUE)$p.value
}

ss14fetmid<-SS2x2(.1,.4,.025,pvalFunc,power=.80,n1start=ceiling(ss14$n.asin),increaseby = 1)





tab<-matrix(NA,3,6,dimnames=list(NULL,
      c("p0","p1","diff","asin","FET","FET mid-p")))
tab[1,]<-c(.3,.6,ss36$n.diff,ss36$n.asin,ss36fet$n1,ss36fetmid$n1)
tab[2,]<-c(.1,.4,ss14$n.diff,ss14$n.asin,ss14fet$n1,ss14fetmid$n1)
tab[3,]<-c(.01,.15,ss15$n.diff,ss15$n.asin,ss15fet$n1,ss15fetmid$n1)

library(xtable)
xtable(ceiling(t(tab)[3:6,]))