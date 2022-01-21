# Section 9.6, simulations on efficiency
library(bpcp)
library(MASS)

nsim<-10^5

N<-30:40

set.seed(034001)
pow<-function(p){ length(p[p<=0.05])/length(p) }

POtest<-function(y1,y2){
  y<-factor(c(y1,y2))
  g<-c(rep(0,length(y1)),rep(1,length(y2)))
  pout<-polr(y~g, contrasts = c(g))
  summary(pout)
}


Pwt<-Pst<-Pwmw<-Ppo<-Pmeld<-rep(NA,length(N))
for (j in 1:length(N)){
  pwt<-pst<-pwmw<-ppo<-pmeld<-rep(NA,nsim)  
  for (i in 1:nsim){
      y1<-rnorm(N[j],sd=9)
      y2<- 7 + rnorm(N[j],sd=9)
      pwt[i]<-t.test(y1,y2)$p.value
      pst[i]<-t.test(y1,y2,var.equal = TRUE)$p.value
      pwmw[i]<-wilcox.test(y1,y2, exact=FALSE)$p.value
      #pmeld[i]<-mdiffmedian.test(y1,y2)$p.value
      #ppo[i]<- POtest(y1,y2)$p.value
  }
  Pwt[j]<-pow(pwt)
  Pst[j]<-pow(pst)
  Pwmw[j]<-pow(pwmw)
  #Pmeld[j]<-pow(pmeld)
}

#d<-data.frame(N,Pst,Pwt,Pwmw,Pmeld)
d<-data.frame(N,Pst,Pwt,Pwmw)
d

EN<-rep(NA,3)
names(EN)<-c("ENst","ENwt","ENwmw")
EN[1]<-approx(Pst,N,0.90)$y
EN[2]<-approx(Pwt,N,0.90)$y
EN[3]<-approx(Pwmw,N,0.90)$y
EN
EN["ENst"]/EN["ENwmw"]
EN["ENst"]/EN["ENwt"]
#library(xtable)
#print(xtable(d,digits=c(0,0,2,3,2,2)),include.rownames=FALSE)
