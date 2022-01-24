# figs 16.4 and 16.5
distributionForms<-function(t,F, excludeTime0=TRUE){
  n<-length(t)
  if (length(F)!=n) stop("length of F must equal length of t")
  if (t[1]!=0) stop("t[1] must be zero")
  delta<- diff(t)
  f<- c(0,diff(F)/delta)
  S<- 1-F
  Sminus<- c(1,S[-n])
  h<- f/Sminus
  H<- cumsum(h*c(0,delta))
  o<- F/S
  #o<-S/F
  D<-data.frame(t=t,F=F,f=f,S=S,Sm=Sminus,h=h,H=H,o=o)
  if (excludeTime0){
    D <- D[-1,]
  }
  D
}





time<- c(0:2000/1000)
time<-c(0:10/10)
d<-distributionForms(time,pexp(time,rate=1))
#d


plotDist<-function(d){
  par(mfrow=c(3,2))
  plot(d$t,d$F,type="l",main="F(t)",lwd=3,xlab="t",ylab="F",ylim=c(0,1))
  plot(d$t,d$f,type="l",main="f(t)",lwd=3,xlab="t",ylab="f")
  plot(d$t,d$S,type="l",main="S(t)",lwd=3,xlab="t",ylab="S",ylim=c(0,1))
  plot(d$t,d$h,type="l",main=expression(lambda(t)),lwd=3,xlab="t",ylab=expression(lambda),ylim=c(0,2*max(d$h)))
  plot(d$t,d$o,type="l",main="Odds(t)=F(t)/S(t-)",lwd=3,xlab="t",ylab="Odds")
  plot(d$t,d$H,type="l",main=expression(Lambda(t)),lwd=3,xlab="t",ylab=expression(Lambda))
  par(mfrow=c(1,1))
}

##plotDist(d[1:2000,])
plotDist(d)


## find MW parameter that goes with proportional odds parameter theta

phi<-function(theta){
  out<-theta*(-1+theta-log(theta))/(1-theta)^2
  out[theta==1]<- 1/2
  out
}





x<-uniroot(rootfunc,c(1,200),MW=.6)


phi(x$root)



Gpo<-function(F,theta){
  F/(F+theta-theta*F)
}

Gph<-function(F,h){
  1-(1-F)^h
}

time<- 0:2000/1000
weirdF<-function(time,frac=1/4){
  F1<-0.2*pexp(time,rate=5) 
  n<-length(time)
  nf<- round(frac*n)
  #F2[1:nf]<- 0 
  #F2[(nf+1):n]<-  0.8*pexp(time[1:(n-nf)],rate=10)
  #F2<-0.8*pexp(time,rate=10)
  F1+F2  
}


weirdF<-function(time){
  set.seed(1)
  n<-length(time)
  y<-rchisq(n,df=1)
  b<-c(round(n/4),round(n/2))
  i1<-1:b[1]
  i2<-(b[1] + 1):b[2]
  i3<-(b[2] + 1):n
  y<-c(sort(y[i1],decreasing=TRUE),sort(y[i2]),sort(y[i3],decreasing=TRUE))
    

  F<- cumsum(y)/sum(y)
  F
}


weirdF2<-function(time){
  n<-length(time)
  b<-c(round(n/4),round(n/2))
  y<-c(b[1]:1,1:(b[2]-b[1]),(b[2]-b[1]):(-n-b[1]+1+2*b[2]))
  y<- y -y[n]+1
  
  F<- cumsum(y)/sum(y)
  F
}
  
  
  
time<-1:2000/1000

#distributionForms(tt,weirdF(tt))

#plot(time,1-weirdF(time))




pophPlots<-function(MW=0.6,dologodds=FALSE,dologhr=FALSE,base="exponential",hazlim=NULL,orlim=NULL,wparms=c(1,3)){
  rootfunc<-function(theta,mw=MW){
    phi(theta) - mw
  }
  OR<- uniroot(rootfunc,c(1e-6,1e6),tol=.Machine$double.eps^0.5)$root
  #OR<-1/OR
  #HR<- MW/(1-MW)
  HR<- (1-MW)/MW
  time<- 0:2000/1000
  if (base=="exponential"){
    F<- pexp(time)
  } else if (base=="log-logistic"){
    F<- plogis(log(time))
  } else if (base=="weirdF"){
    F<- weirdF(time)
  } else if (base=="Weibull"){
    F<-pweibull(time,shape=wparms[1],scale=wparms[2])
  }
  GPO<-Gpo(F,OR)
  GPH<-Gph(F,HR)
  par(mfrow=c(3,2))
  S<-1-F
  TPO<- 1- GPO
  TPH<-1-GPH
  dF<-distributionForms(time,F)
  dGPO<-distributionForms(time,GPO)
  dGPH<-distributionForms(time,GPH)
  
  par(mfrow=c(3,2))
  par(mar=c(4.1,4.1,0,2.1))
  par(oma=c(0,0,4,0))
  
  LWD<-c(5,3)
  COL<-c("black","gray")
  
  plot(time,S,xlab="",ylab="S",type="l",main="",lwd=LWD[1],col=COL[1],ylim=c(0,1))
  lines(time,TPH,lwd=LWD[2],col=COL[2])
 
  plot(time,S,xlab="",ylab="S",type="l",main="",lwd=LWD[1],col=COL[1], ylim=c(0,1))
  lines(time,TPO,lwd=LWD[2],col=COL[2])
  if (dologhr){
    
    # calculate log(hazard) range
    logHlim<-range(log(c(dF$h,dGPH$h,dGPO$h)))
    
    plot(dF$t,log(dF$h),xlab="",ylab="log(hazard)",type="l",main="",lwd=LWD[1],col=COL[1],ylim=logHlim)
    lines(dGPH$t,log(dGPH$h),lwd=LWD[2],col=COL[2])
    
    
    plot(dF$t,log(dF$h),xlab="",ylab="log(hazard)",type="l",main="",lwd=LWD[1],col=COL[1],ylim=logHlim)
    lines(dGPO$t,log(dGPO$h),lwd=LWD[2],col=COL[2])
  } else {
    
    # calculate log(hazard) range
    if (is.null(hazlim)){
      Hlim<-range((c(dF$h,dGPH$h,dGPO$h)))+c(-.1,+.1)      
    } else {
      Hlim<- hazlim
    }

    
    plot(dF$t,(dF$h),xlab="",ylab="hazard",type="l",main="",lwd=LWD[1],col=COL[1],ylim=Hlim)
    lines(dGPH$t,(dGPH$h),lwd=LWD[2],col=COL[2])
    
    
    plot(dF$t,(dF$h),xlab="",ylab="hazard",type="l",main="",lwd=LWD[1],col=COL[1],ylim=Hlim)
    lines(dGPO$t,(dGPO$h),lwd=LWD[2],col=COL[2])
  }
   
  if (dologodds){
    
    # calculate log(odds) range
    logOddsLim<-range(log(1/c(dF$o,dGPH$o,dGPO$o)))
    
    
    plot(dF$t,log(1/dF$o),xlab="time",ylab="-log(odds)",type="l",main="",lwd=LWD[1],col=COL[1],
         ylim=logOddsLim)
    lines(dGPH$t,log(1/dGPH$o),lwd=LWD[2],col=COL[2])
    
    plot(dF$t,log(1/dF$o),xlab="time",ylab="-log(odds)",type="l",main="",lwd=LWD[1],col=COL[1],
         ylim=logOddsLim)
    lines(dGPO$t,log(1/dGPO$o),lwd=LWD[2],col=COL[2])
  } else {
    if (is.null(orlim)){
      # calculate (odds) range
      OddsLim<-range(c(dF$o,dGPH$o,dGPO$o)) 
    } else {
      OddsLim<-orlim
    }
   
    
    plot(dF$t,(dF$o),xlab="time",ylab="odds",type="l",main="",lwd=LWD[1],col=COL[1],
         ylim=OddsLim)
    lines(dGPH$t,(dGPH$o),lwd=LWD[2],col=COL[2])
    
    plot(dF$t,(dF$o),xlab="time",ylab="odds",type="l",main="",lwd=LWD[1],col=COL[1],
         ylim=OddsLim)
    lines(dGPO$t,(dGPO$o),lwd=LWD[2],col=COL[2])
  }
   
  mtext("Proportional Hazards",side=3,line=1,adj=0.2,outer=TRUE)
  
  mtext("Proportional Odds",side=3,line=1,adj=.8,outer=TRUE)
  
  
  par(mfrow=c(1,1),mar=c(5,4,4,2)+.1,oma=c(0,0,0,0)) 
}

# Figure 16.4
pophPlots(.6, base="exponential")
#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenPlotPOPHexpoential.pdf")


#pophPlots(.6,base="log-logistic")
#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenPlotPOPHloglogistic.pdf")

# Figure 16.5
pophPlots(.6,base="weirdF",hazlim=c(0,20),orlim=c(0,20))
#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenPlotPOPHweird.pdf")

#pophPlots(.6,base="Weibull",wparms=c(3,1),orlim=c(0,12))
#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenPlotPOPHweibull31.pdf")

