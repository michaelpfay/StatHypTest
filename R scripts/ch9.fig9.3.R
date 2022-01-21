# Figure 9.3

# See ch9.example9.4.WMWMixBiNorm.R for how I found the parameter estimates 
# need to source it from some functions also
source("ch9.example9.4.R")
library(mvtnorm)

plotMixBiNorm<-function(ug1,ug2,vg1,vg2,rhog,
                        uh1,uh2,vh1,vh2,rhoh, lambda,ngrid=500,...){
    # get range
    ff<-3
    min1<- min(ug1-ff*sqrt(vg1), uh1-ff*sqrt(vh1))
    min2<- min(ug2-ff*sqrt(vg2), uh2-ff*sqrt(vh2))
    max1<- max(ug1+ff*sqrt(vg1), uh1+ff*sqrt(vh1))
    max2<- max(ug2+ff*sqrt(vg2), uh2+ff*sqrt(vh2))

    MIN<-min(min1,min2)
    MAX<-max(max1,max2)
    
    #x<-seq(min1,max1,length.out=ngrid)
    #y<-seq(min2,max2,length.out=ngrid)
    x<-seq(MIN,MAX,length.out=ngrid)
    y<-seq(MIN,MAX,length.out=ngrid)
    z<-matrix(NA,ngrid,ngrid)
    for (i in 1:ngrid){
        for (j in 1:ngrid){
            z[i,j]<- lambda*dmvnorm(c(x[i],y[j]),
                   mean=c(ug1,ug2),
                   sigma=matrix(c(vg1,rep(rhog*sqrt(vg1*vg2),2),vg2),2,2)) +
                 (1-lambda)*dmvnorm(c(x[i],y[j]),
                   mean=c(uh1,uh2),
                   sigma=matrix(c(vh1,rep(rhoh*sqrt(vh1*vh2),2),vh2),2,2))
        }
    }
    #persp(x,y,z,...)
    return(list(x=x,y=y,z=z))
}


#out<-plotMixBiNorm(ug1=2,ug2=3,vg1=1,vg2=1,rhog=.9,
#              uh1=5,uh2=0,vh1=1,vh2=.1,rhoh=0,lambda=.65,ngrid=100)
out<-plotMixBiNorm(ug1=3,ug2=2,vg1=1,vg2=1,rhog=.9,
                   uh1=0,uh2=5,vh1=.1,vh2=1,rhoh=0,lambda=.65,ngrid=100)

par(mfrow=c(1,2))

contour(out$x,out$y,out$z,ylim=c(-1,7),xlim=c(-2,6.5),levels=c(.001,0.01,.05,.1,.15,.2),xlab="Y(1)",ylab="Y(2)")
lines(c(-10,10),c(-10,10),lty=2,lwd=3,col=gray(.5))

plotMixBiNormMarginals<-function(ug1,ug2,vg1,vg2,rhog,
                        uh1,uh2,vh1,vh2,rhoh, lambda,ngrid=100,...){
  # get range
  ff<-3
  min1<- min(ug1-ff*sqrt(vg1), uh1-ff*sqrt(vh1))
  min2<- min(ug2-ff*sqrt(vg2), uh2-ff*sqrt(vh2))
  max1<- max(ug1+ff*sqrt(vg1), uh1+ff*sqrt(vh1))
  max2<- max(ug2+ff*sqrt(vg2), uh2+ff*sqrt(vh2))
  
  MIN<-min(min1,min2)
  MAX<-max(max1,max2)
  
  #x<-seq(min1,max1,length.out=ngrid)
  #y<-seq(min2,max2,length.out=ngrid)
  x<-y<-seq(MIN,MAX,length.out=ngrid)
  fx<- lambda*dnorm(x,mean=ug1,sd=sqrt(vg1)) + (1-lambda)*dnorm(x,uh1,sd=sqrt(vh1))
  fy<- lambda*dnorm(y,mean=ug2,sd=sqrt(vg2)) + (1-lambda)*dnorm(y,uh2,sd=sqrt(vh2))
  plot(x,fx,type="l",lty=1,lwd=4,col=gray(.5),...)
  lines(y,fy,lty=2,lwd=3)
  legend("topright",legend=c("Marginal for Y(1)","Marginal for Y(2)"),lty=c(1,2),lwd=c(4,3),col=c(gray(.5),"black"))
  return(list(x=x,fx=fx,y=y,fy=fy))
  
}

#plotMixBiNormMarginals(ug1=2,ug2=3,vg1=1,vg2=1,rhog=.9,
#                       uh1=5,uh2=0,vh1=1,vh2=.1,rhoh=0,lambda=.65,ylim=c(0,.5),xlab="y",ylab="f(y)")

plotMixBiNormMarginals(ug1=3,ug2=2,vg1=1,vg2=1,rhog=.9,
                       uh1=0,uh2=5,vh1=.1,vh2=1,rhoh=0,lambda=.65,ylim=c(0,.5),xlab="y",ylab="f(y)")

MixBiNorm(c(3,2),1,1,.9,c(0,5),.1,1,0,0.65)

#dev.print(pdf,file="../graph/ch2OrdMixBiNorm.pdf")
