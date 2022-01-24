# Fig 19.1
# plot simulation from 3 distributions
library(VGAM)
set.seed(40201)
n<-100
x.norm<- 5 + 3* rnorm(n)
x.t<- -3 + 0.2*rt(n,df=5)
x.lgamma<- rlgamma(n,shape=0.2,location=124,scale=30)

par(mfrow=c(2,3))
hist(x.norm,main="Normal")
hist(x.t,main="t df=5")
hist(x.lgamma,main="log-gamma")

qqnorm(x.norm)
qqline(x.norm)


qqnorm(x.t)
qqline(x.t)

qqnorm(x.lgamma)
qqline(x.lgamma)

#dev.print(pdf,file="../book/graph/chNI.qqnormPlots.pdf")
