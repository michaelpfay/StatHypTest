# Fig 6.5
library(mvtnorm)
set.seed(10301)
sim<-function(rho,nsim=10^5){
    x<-rmvnorm(nsim,sigma=matrix(c(1,rho,rho,1),2,2))
    cor(pnorm(x[,1]),pnorm(x[,2]))
}    


RHO<- -100:100/100
rhat<-rep(NA,length(RHO))
for (i in 1:length(RHO)){
    rhat[i]<-sim(RHO[i])
}

plot(RHO,rhat,type="n",xlab="Pearson Correlation",ylab="Association")
lines(c(-1,1),c(-1,1),lwd=3,lty=2,col="gray")
s<-supsmu(RHO,rhat)
lines(RHO,s$y,lwd=3,col="black")

lines(RHO,(2/pi)*asin(RHO),lwd=3,lty=3)
#legend("bottomright",legend=c("Kendall's tau","limit E(Spearman Correlation)","Pearson Correlation (identity)"),
#     lwd=c(3,3,3),lty=c(3,1,2),col=c("black","black","gray"))

#dev.print(pdf,"../graph/chPaiNormalCorrelations.pdf")
