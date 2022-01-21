# Section 13.1 (unnumbered table)


# plot FWER and FDR for m tests that are correlated with rho=r
par(mfrow=c(1,1))
library(mvtnorm)
alpha<-0.05


pmaxZ<-function(x,m,rho){
  if (m==1){
    out<-1- integrate(dnorm,-x,x)$value
  } else {
    R<-matrix(rho,m,m)+ diag(rep(1-rho,m))
    out<-1-pmvnorm(lower=rep(-x,m),upper=rep(x,m),corr=R)    
  }
  out
}


m<-50
FWER2<-rep(NA,m)

rho<-.5
M<-c(2,3,5,10,20,100)
m<-length(M)
FWER2<-rep(NA,m)
for (i in 1:m){
  FWER2[i]<-pmaxZ(qnorm(1-alpha/2),M[i],rho)
}

FWER<-1-(1-alpha)^(M)

round(cbind(FWER,FWER2),4)

plot(1:m,FWER,type="l",ylim=c(0,1),xlab="m",ylab="FWER")
lines(1:m,FWER2,lty=2,lwd=3,col=gray(.8))

## FDR, it is easy with rho=0
# Dj ~ Bernoulli(alpha)
# sum Dj =X ~ Binomial(m, alpha)
sum( c(0:m) * dbinom(0:m,m,alpha) )



