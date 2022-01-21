# Figure 10.1


# plot size from Wald, Score and LRT for the one-sample normal test of mean with variance unspecified
# for derivation see p. 134 of Boos and Stefanski
n<- 5:100
alpha<-0.05
Chi<- qchisq(1-alpha,1)
sizeWald<- 1- pf( ((n-1)/n)*Chi,1,n-1)
sizeWald
sizeScore<- 1- pf( ((n-1)/(n-Chi))*Chi,1,n-1)
sizeScore

sizeLRT<- 1- pf( (n-1)*(exp(Chi/n)-1),1,n-1)
sizeLRT




par(mfrow=c(1,2))

plot(n,sizeWald,type="l",ylim=c(0,.1),lwd=2,ylab="size",main="One-sample Normal Problem")
lines(n,sizeScore,lty=2,lwd=3)
lines(n,sizeLRT,lty=3,lwd=3)
lines(n,rep(alpha,length(n)),lwd=3,col=gray(.7))

legend(60,.1,legend=c("Wald","LRT","score"),lwd=c(2,3,3),lty=c(1,3,2))



### Simulation for logistic regression
### Erica did the simulation with SAS
### see /sas/chGMScoreLogisticSim.sas
n2<- c(10,15,25,50)
sizeWald2<- c(0,.006,.03,.042)
sizeScore2<-c(0.056,.053,.051,.051)
sizeLRT2<-c(0.095,.073,.062,.055)

plot(n2,sizeWald2,type="l",ylim=c(0,.1),lwd=2,ylab="size",main="One-sample Logistic Problem")
lines(n2,sizeScore2,lty=2,lwd=3)
lines(n2,sizeLRT2,lty=3,lwd=3)
lines(n2,rep(alpha,length(n2)),lwd=3,col=gray(.7))

legend(32,.1,legend=c("Wald","LRT","score"),lwd=c(2,3,3),lty=c(1,3,2))




#dev.print(pdf,file="../graph/chGMAsymOneSampleNormal2.pdf")