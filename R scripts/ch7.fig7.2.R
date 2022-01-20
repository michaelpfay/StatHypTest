# Figure 7.2
# programs may take several minutes to run
# check type I error rates 
library(exact2x2)


plotPowersFisher<-function(N1,N2, ngrid=100){
  theta<- 0:ngrid/ngrid
  
  fet<-fmidp<-fmidpmeld<-rep(NA,length(theta))
  
  for (i in 1:length(theta)){
    fmidp[i]<-Power2x2(N1,N2,theta[i],theta[i],0.025, pvalFunc=
                         function(x1,n1,x2,n2){
                           exact2x2(matrix(c(x1,n1-x1,x2,n2-x2),2,2), alternative="greater", 
                                    midp=TRUE)$p.value
                         }
    )
    fet[i]<-Power2x2(N1,N2,theta[i],theta[i],0.025, pvalFunc=
                       function(x1,n1,x2,n2){
                         fisher.test(matrix(c(x1,n1-x1,x2,n2-x2),2,2), alternative="greater")$p.value
                       }
    )
    fmidpmeld[i]<-Power2x2(N1,N2,theta[i],theta[i],0.025, pvalFunc=
                       function(x1,n1,x2,n2){
                         binomMeld.test(x1,n1,x2,n2, conf.int=FALSE, midp=TRUE, alternative="greater")$p.value
                       }
    )
  }
  
  
  
  #cbind(fmidp,fet)
  
  
  
  
  #plot(theta,fmidp,type="l",ylim=c(0,0.035),main=paste("n1=",N1," n2=",N2),xlab=expression(theta),ylab="Pr[ reject ]")
  plot(theta,fmidp,type="l",ylim=c(0,0.035),main=bquote(n[1] == .(N1)~~~n[2]==.(N2)),
       xlab=expression(theta),ylab="Pr[ reject ]")
  lines(c(-1,2),c(0.025,0.025),lty=2)
  lines(theta,fet,lwd=4,col=gray(.7))
  lines(theta,fmidpmeld,lwd=4,lty=2,col=gray(.7))
  lines(theta,fmidp,lwd=1,lty=1,col="black")
}







par(mfrow=c(2,2))

plotPowersFisher(10,10)
plotPowersFisher(7,12)
plotPowersFisher(20,20)
plotPowersFisher(35,25)

#dev.print(pdf,file="../graph/ch2binMidpFisherVSCentralFisher.pdf")
#dev.print(pdf,file="ch2binMidpFisherVSCentralFisher.pdf")