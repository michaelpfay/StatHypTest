# Figs 17.2 and 17.3 and 17.4
PI<- function(w){
  plogis(6*w)
}
W<- c(-20:20/10)

PI2<- function(w){
  pnorm(3*w)
}

par(mfrow=c(2,1))
plot(W,W, type="l",lwd=2,ylab=expression(paste("E( Y | ","W )")))
m2<-W
m2[W<0]<-0
lines(W,m2,col="gray",lwd=6,lty=2)
lines(W,W,lwd=2)
legend("topleft",lwd=c(2,6),col=c("black","gray"),lty=c(1,2),legend=c("Scenario 1","Scenario 2"))



plot(W,PI(W), type="l",lwd=2,ylab=expression(pi(W)))
lines(W,PI2(W),col="gray",lwd=6,lty=2)
lines(W,PI(W),lwd=2)
legend("topleft",lwd=c(2,6),col=c("black","gray"),lty=c(1,2),legend=c("Scenario 1","Scenario 2"))

#dev.print(pdf, file="../book/graph/chMiNormalSim1.pdf")


##########################################
# Create and Plot 2 simulated data sets
##########################################

N<-100
SIGW<-1
SIGE<-1

dorow<-function(Y,W,R){
  lout<-lm(Y~W,subset=(R==1))
  yhat<- predict(lout,newdata=data.frame(W=W))
  muobs<- mean(Y[R==1])
  muRI<- mean(yhat)
  # IPW using logistic regression
  gout<-glm(R~W,family=binomial)
  pihat<-predict(gout,type="response")
  # Order the W so that we can plot the line
  muIPW<-mean(R*Y/pihat)
  muIPWb<- sum(R*Y/pihat)/sum(R/pihat)
  piBeta<- coef(gout)
  
  out<-list(mu=c(ALL=mean(Y),OBS=muobs,
            RI=muRI,
            IPW=muIPW),piBeta=piBeta)
  out
}

createData<-function(scenario,n=N){
  W<-rnorm(n,sd=SIGW)
  delta<-function(x){ y<-x; y[x<0]<-0; y }
  if (scenario==1){
    Y<- W + rnorm(n,sd=SIGE)    
  } else {
    Y<- delta(W) + rnorm(n,sd=SIGE)
  }

  if (scenario==1){
    p<-PI(W)
  } else {
    p<-PI2(W)
  }
  R<-rbinom(n,1,p)
  out<-list(Y=Y,W=W,R=R)
  out
}

EYgW<-function(scenario,sigw=SIGW,sige=SIGE,nMC=1e5){
  if (scenario==1){
    out<-0
  } else {
    delta<-function(x){ y<-x; y[x<0]<-0; y }
    out<- mean( delta(rnorm(nMC,sd=sigw))+ rnorm(nMC,sd=sige))
  }
  out
}

EYgW2<-EYgW(2,nMC=1e7)
EYgW2


set.seed(40921)
d1<-createData(1)
set.seed(3901988)
d2<-createData(2)

par(mfrow=c(3,2))
rangeW<-range(c(d1$W,d2$W))
rangeY<-range(c(d1$Y,d2$Y))

plot(d1$W,d1$Y,type="n",xlim=rangeW,ylim=rangeY,xlab="W",ylab="Y",main="Scenario 1")
points(d1$W[d1$R==1],d1$Y[d1$R==1])
plot(d2$W,d2$Y,type="n",xlim=rangeW,ylim=rangeY,xlab="W",ylab="Y",main="Scenario 2")
points(d2$W[d2$R==1],d2$Y[d2$R==1])

hist(d1$W[d1$R==1],xlim=rangeW,xlab="W when R=1",main="Histogram of W for R=1")
hist(d2$W[d2$R==1],xlim=rangeW,xlab="W when R=1",main="Histogram of W for R=1")


hist(d1$W[d1$R==0],xlim=rangeW,xlab="W when R=0",main="Histogram of W for R=0")
hist(d2$W[d2$R==0],xlim=rangeW,xlab="W when R=0",main="Histogram of W for R=0")

#dev.print(pdf, file="../book/graph/chMiNormalSim2.pdf")

nsim<-1e3

COEF1<-COEF2<- matrix(NA,nsim,2)
muRes1<-muRes2<-matrix(NA,nsim,4)



for (i in 1:nsim){
  d1<-createData(1)
  d2<-createData(2)
  x1<-dorow(d1$Y,d1$W,d1$R)
  x2<-dorow(d2$Y,d2$W,d2$R)
  COEF1[i,]<- x1$piBeta
  muRes1[i,]<-x1$mu
  COEF2[i,]<- x2$piBeta
  muRes2[i,]<-x2$mu
}

dimnames(muRes1)[[2]]<-dimnames(muRes2)[[2]]<- names(x1$mu)


#hist(muResults[,1],xlim=c(-10,10),main=dimnames(muResults)[[2]][1])
#hist(muResults[,2],xlim=c(-10,10),main=dimnames(muResults)[[2]][2])
#hist(muResults[,3],xlim=c(-10,10),main=dimnames(muResults)[[2]][3],breaks=c(-Inf,-20:20/2,Inf))
#hist(muResults[,4],xlim=c(-10,10),main=dimnames(muResults)[[2]][4],breaks=c(-Inf,-20:20/2,Inf))


library(reshape2)
library(ggplot2)
library(gridExtra)
md1<-melt(muRes1,value.name="estimate")
md2<-melt(muRes2,value.name="estimate")
#par(mfrow=c(1,1))
#plot(mdata$estimate~mdata$Var2)
#lines(c(-10,10),c(0,0),col="gray",lwd=3,lty=2)



# Basic violin plot
ysimRange<-range(c(range(muRes1),range(muRes2)))
p1 <- ggplot(md1, aes(x=Var2, y=estimate)) + 
  geom_violin() + scale_y_continuous(limits=ysimRange) +
  geom_hline(yintercept = 0)+ggtitle("Scenario 1")+xlab("")
#+  geom_jitter(shape=16, position=position_jitter(0.2))
#+ geom_dotplot(binaxis='y', stackdir='center', dotsize=.01)
p2 <- ggplot(md2, aes(x=Var2, y=estimate)) + 
  geom_violin() + scale_y_continuous(limits=ysimRange)+
  geom_hline(yintercept = EYgW2)+ggtitle("Scenario 2")+xlab("")

grid.arrange(p1,p2,nrow=1)

#dev.print(pdf, file="../book/graph/chMiNormalSim3.pdf")
