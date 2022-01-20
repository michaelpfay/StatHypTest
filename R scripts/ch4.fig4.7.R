library(exactci)

powerBinom<-function(n,theta1=0.7,theta0=.5,alpha=0.025,one.sided=TRUE,...){
    reject<-rep(FALSE,n+1)
    if (one.sided){
       for (i in 1:(n+1)){
          alt<-ifelse(theta1>theta0,"greater","less")
          if (binom.test(i-1,n,p=theta0,alternative=alt)$p.value<=alpha) reject[i]<-TRUE
       }
    } else {
       for (i in 1:(n+1)){
          if (binom.exact(i-1,n,p=theta0,...)$p.value<=alpha) reject[i]<-TRUE
       }
    }
    list(power=sum(dbinom(0:n,n,theta1)[reject]),reject=reject)
}

N<-2:100
power1<-power2<-power2b<-rep(NA,length(N))
for (i in 1:length(N)){
    power1[i]<-powerBinom(N[i],alpha=0.025,one.sided=TRUE)$power
    power2[i]<-powerBinom(N[i],alpha=0.05,one.sided=FALSE,tsmethod="central")$power
    power2b[i]<-powerBinom(N[i],alpha=0.05,one.sided=FALSE,tsmethod="blaker")$power
}
names(power1)<-names(power2)<-N
p61<-powerBinom(61,alpha=)
p62<-powerBinom(62)

range(power2b-power2)


# Figure 4.7
plot(N,power1,type="l",lwd=6,col=gray(.8),ylab="Power")
lines(N,power2,lwd=2,col="black")
lines(c(-10,200),c(.8,.8),lty=2,col=gray(.5))
lines(c(49,49),c(power1[N=49],-10),lty=2,col=gray(.5))
lines(c(54,54),c(power1[N=54],-10),lty=2,col=gray(.5))
axis(1,at=c(49,54),labels=c(49,54),cex.axis=.8)

dev.print(pdf,file="../graph/ch1bPowerBinom.pdf")
