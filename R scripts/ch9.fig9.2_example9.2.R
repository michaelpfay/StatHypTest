# Fig 9.2 and Example 9.2
library(perm)
days<-1:1200
a1<-50
b1<-1000
a2a<- 50
b2a<-1100
a2b<-1
b2b<-30
pout<-.1



par(mfrow=c(1,1))

f1<-dweibull(days,a1,b1)
f2a<-dweibull(days,a2a,b2a)
f2b<- dweibull(days,a2b,b2b)
f2<- pout*f2b + (1-pout)*f2a
plot(c(days,days),c(f1,f2),type="n",xlab="days",ylab="",axes=FALSE)
axis(1)
box()
lines(days,f2,lwd=6,col=gray(.8))
lines(days,f1,lwd=2)


## mean of f1
mu1<- b1*gamma(1+1/a1)
mu1
qweibull(c(.25,.75),a1,b1)

## mean of f2
mu2<- (1-pout)* b2a*gamma(1+1/a2a)+pout*b2b*gamma(1+1/a2b)
mu2a<- b2a*gamma(1+1/a2a)
mu2b<-b2b*gamma(1+1/a2b)
mu2a
mu2b
mu2
qweibull(c(.25,.75),a2a,b2a)
qweibull(c(.25,.75),a2b,b2b)



nsim<-10^4
ni<-5
chm<-chooseMatrix(2*ni,ni)
dim(chm)
pval<-rep(NA,nsim)
for (i in 1:nsim){
    y1<-rweibull(ni,a1,b1)
    y2<-rweibull(ni,a2a,b2a)
    y2b<-rweibull(ni,a2b,b2b)
    pick<- rbinom(ni,1,pout)==1
    if (any(pick)) y2[pick]<- y2b[pick]
    pval[i]<-permTS(y1,y2,exact=true,method="exact.ce",alternative="less",control=permControl(cm=chm))$p.value
}

length(pval[pval<=0.05])/nsim
dbinom(0,ni,pout)


#dev.print(pdf,file="../graph/ch2OrdWeibull.pdf")

### do a simulated data example
set.seed(1994839)
y1<-rweibull(ni,a1,b1)
y2<-rweibull(ni,a2a,b2a)
y1<-round(y1,0)
y2<-round(y2,0)
y1
y2
calcCI<-function(y1,y2,D){
    pL<-pU<-rep(NA,length(D))
    for (i in 1:length(D)){
        pL[i]<-permTS(y2-D[i],y1,exact=true,method="exact.ce",
             alternative="greater",control=permControl(cm=chm))$p.value
        pU[i]<-permTS(y2-D[i],y1,exact=true,method="exact.ce",
             alternative="less",control=permControl(cm=chm))$p.value
    }
    CRegion<-D[pL>0.025 & pU>0.025]
    CI<-c(min(CRegion),max(CRegion))
    c(pL0=pL[D==0],pU0=pU[D==0],lower=CI[1],upper=CI[2])
}

D<- -100:200
result<-calcCI(y1,y2,D)
result
D<- c(0,760:2000)/10
result<-calcCI(y1,y2,D)
result
D<-c(0,7690:7710,15200:15300)/100
result<-calcCI(y1,y2,D)
result
1/252

log10(y1)
log10(y2)
D<- c(0,3200:3300,6300:6500)/100000
result<-calcCI(log10(y1),log10(y2),D)
result
10^result[3:4]


