# Section 12.2.2 simulations
library(lme4)

t0<-proc.time()

createBetaBinomial<-function(m,a,b){
  n<-length(m)
  pi<-rbeta(n,a,b)
  x<-rbinom(n,m,pi)
  list(x=x,m=m)
}

nsim<-10000
#nsim<-100
N<- 50
A1<-1
B1<-3
A2<-1
B2<-3
#A2<- 1.2
#B2<- 2.8
### Simulation on Type I error rate
###    Sample size random and not related to response
set.seed(101)
options(warn=0)

pvalQB<-pvalB<-pvalt<-phi<-pvalmix<-rep(NA,nsim)
for (i in 1:nsim){
  #print(paste("i=",i))
  out1<-createBetaBinomial(1+rpois(N,10),A1,B1)
  out2<-createBetaBinomial(1+rpois(N,10),A2,B2)
  
  x<-c(out1$x,out2$x)
  m<-c(out1$m,out2$m)
  y<- x/m
  g<-c(rep(0,N),rep(1,N))
  
  pvalt[i]<-t.test(y~g)$p.value
  
  gqb<-glm(cbind(x,m-x)~g,family=quasibinomial())
  pvalQB[i]<-summary(gqb)$coef["g","Pr(>|t|)"]
  phi[i]<-summary(gqb)$dispersion
  
  gb<-glm(cbind(x,m-x)~g,family=binomial())
  pvalB[i]<-summary(gb)$coef["g","Pr(>|z|)"]
  
  id<-1:length(x)
  gmix<-try(glmer(cbind(x,m-x)~ g+(1|id), family=binomial))
  if (class(gmix)!="try-error"){
    if (length(summary(gmix)$optinfo$conv$lme4)==0){
      pvalmix[i]<-summary(gmix)$coef["g","Pr(>|z|)"]    
    }
  }
}

length(phi[phi<1])/nsim
length(pvalt[pvalt<=0.05])/nsim
length(pvalB[pvalB<=0.05])/nsim
length(pvalQB[pvalQB<=0.05])/nsim
any(is.na(pvalmix))
length(pvalmix[is.na(pvalmix)])/nsim
# set 
Imiss<- is.na(pvalmix)
nMiss<- length(pvalmix[Imiss])
pvalmixNoMiss<-pvalmixFix<-pvalmix
pvalmixFix[Imiss]<-1
pvalmixNoMiss<- pvalmix[!is.na(pvalmix)]
length(pvalmixFix[pvalmixFix<=0.05])/nsim
length(pvalmixNoMiss[pvalmixNoMiss<=0.05])/length(pvalmixNoMiss)

t1<-proc.time()
t1-t0

#####################################################
## Simulation where sample size is related to 
## response
##
####################################################

createSSrelResp<-function(n,a,b){
  
  pi<-rbeta(n,a,b)
  #m<- round(10^(1+pi))
  m<- 1+rpois(n,10^(1+pi))
  x<-rbinom(n,m,pi)
  list(x=x,m=m)
}

set.seed(553)

pvalQB<-pvalB<-pvalt<-phi<-pvalmix<-rep(NA,nsim)
for (i in 1:nsim){
  out1<-createSSrelResp(N,A1,B1)
  out2<-createSSrelResp(N,A2,B2)
  
  x<-c(out1$x,out2$x)
  m<-c(out1$m,out2$m)
  y<-x/m
  g<-c(rep(0,N),rep(1,N))
  
  pvalt[i]<-t.test(y~g)$p.value
  
  
  gqb<-glm(cbind(x,m-x)~g,family=quasibinomial())
  pvalQB[i]<-summary(gqb)$coef["g","Pr(>|t|)"]
  phi[i]<-summary(gqb)$dispersion
  
  gb<-glm(cbind(x,m-x)~g,family=binomial())
  pvalB[i]<-summary(gb)$coef["g","Pr(>|z|)"]
  id<-1:length(x)
  gmix<-try(glmer(cbind(x,m-x)~ g+(1|id), family=binomial))
  if (class(gmix)!="try-error"){
    if (length(summary(gmix)$optinfo$conv$lme4)==0){
      pvalmix[i]<-summary(gmix)$coef["g","Pr(>|z|)"]    
    }
  }
}

length(phi[phi<1])/nsim
length(pvalt[pvalt<=0.05])/nsim
length(pvalB[pvalB<=0.05])/nsim
length(pvalQB[pvalQB<=0.05])/nsim
any(is.na(pvalmix))
length(pvalmix[is.na(pvalmix)])/nsim
# set 
Imiss<- is.na(pvalmix)
nMiss<- length(pvalmix[Imiss])
pvalmixNoMiss<-pvalmixFix<-pvalmix
pvalmixFix[Imiss]<-1
pvalmixNoMiss<- pvalmix[!is.na(pvalmix)]
length(pvalmixFix[pvalmixFix<=0.05])/nsim
length(pvalmixNoMiss[pvalmixNoMiss<=0.05])/length(pvalmixNoMiss)


