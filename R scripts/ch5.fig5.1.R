# Figure 5.1
source("ch5.ttestSimGetLogGammaParms.R")

## plot simulated distributions
## standardize so that F(.25)=-1, and F(.75)=1


plotDist<-function(
    x= (-500:500)/100,
    tick=c(-1,0,1),
    dfunc=dnorm,
    TITLE="Normal"){

    y<-dfunc(x)
    ytick=dfunc(tick)
    plot(x,y,type="l",xlab="",ylab="",main=TITLE,lwd=3,yaxs="i",ylim=c(0,.4),
                 xlim=c(-5,5),xaxs="i")
    segments(tick[c(1,3)],c(-1,-1),tick[c(1,3)],ytick[c(1,3)],lty=3,lwd=3,col=gray(.5))
    segments(tick[c(2)],c(-1),tick[c(2)],ytick[c(2)],lty=1,lwd=3,col=gray(.5))

    lines(x,y,lwd=3)
}


par(mfrow=c(2,2),mar=c(2.5,3,3,1))
## Normal
dNormal<-function(x){ dnorm(x,sd=1/qnorm(.75)) }
##check
#pnorm(c(-1,0,1),sd=1/qnorm(.75))
plotDist(dfunc=dNormal,TITLE="Normal")
## t
dtdist<- function(x,df=1){ 
    q<-qt(.75,df)
    q*dt(q*x,df)
}
## check
#integrate(dtdist,-Inf,-1)
#integrate(dtdist,-1,1)
plotDist(dfunc=dtdist,TITLE="Cauchy")

## t ddf=2
dist<-function(x){ dtdist(x,df=2) }
plotDist(dfunc=dist,TITLE="t df=2")

## log-gamma


plotDist(dfunc=logGammaDist,tick=c(-1,logGammaMean,1),
   TITLE=paste0("log gamma (shape=",SHAPE,")"))

#dev.print(pdf,"../graph/ch1OrdTtestSimDistns.pdf")

