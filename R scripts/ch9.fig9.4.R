# Figure 9.4 + simulation on efficiency
# Made up data
n<-112




par(mfrow=c(1,2))
set.seed(1)
y<- 10*rt(n,df=17) + 1000 
mean(y)
sd(y)
hist(y, xlab="lifetime (in days)", main="(a)")

set.seed(1)
y<-rweibull(n,shape=50, scale=1000)
mean(y)
sd(y)
hist(y, xlab="lifetime (in days)", main="(b)")


#dev.print(pdf,file="../graph/ch2OrdMiceHistogram.pdf")

power.t.test(delta=7,sd=9, sig.level=0.025, power=.9, type="two.sample",
             alternative="one.sided")

nsim<-10^5
set.seed(1010)
N<-c(30:40)
std<-9
shift<-7
getpow<-function(p){ length(p[p<=0.025])/length(p)}
powwmw<-powt<-powtw<-rep(NA,length(N))
for (h in 1:length(N)){
  pwmw<-pt<-ptw<-rep(NA,nsim)
  for (i in 1:nsim){
     pwmw[i]<-wilcox.test(std*rnorm(N[h]),std*rnorm(N[h])+shift,alternative="less")$p.value
     pt[i]<-t.test(std*rnorm(N[h]),std*rnorm(N[h])+shift,alternative="less",var.equal = TRUE)$p.value
     ptw[i]<-t.test(std*rnorm(N[h]),std*rnorm(N[h])+shift,alternative="less", var.equal=FALSE)$p.value
  }
     powwmw[h]<- getpow(pwmw)
     powt[h]<-getpow(pt)
     powtw[h]<-getpow(ptw)  
  }

names(powwmw)<-names(powt)<-names(powtw)<-paste("N1=",N)

powwmw
powt
powtw

approx(powwmw,N,xout=.9)$y
approx(powt,N,xout=.9)$y
approx(powtw,N,xout=.9)$y