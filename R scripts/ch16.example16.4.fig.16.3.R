# example 16.4, fig 16.3
library(survival)
library(bpcp)
#library(ggplot2)
data(leuk2)
levels(leuk2$treatment)
leuk2$treatment<-factor(leuk2$treatment,levels=c("placebo","6-MP"))
sr<-survreg(Surv(time,status)~treatment, data=leuk2)

pct<-1:9999/10000
ptime<-predict(sr,newdata=data.frame(treatment="placebo"),type="quantile",p=pct)
ttime<-predict(sr,newdata=data.frame(treatment="6-MP"),type="quantile",p=pct)

km<-survfit(Surv(time,status)~treatment, data=leuk2)



plot(km,lwd=2,lty=c(1,2))


# Add cox regression lines
levels(leuk2$treatment)
leuk2$treatment<-factor(leuk2$treatment,levels=c("placebo","6-MP"))
cph<-coxph(Surv(time,status)~treatment, data=leuk2)
cfit<-survfit(cph,newdata=data.frame(treatment=factor(levels(leuk2$treatment))))
lines(cfit,col="gray",lty=c(1,2),lwd=6)
# Weibull lines
lines(ptime,1-pct,col="gray",lwd=2,lty=1)
lines(ttime,1-pct,col="gray",lwd=2,lty=2)
# put KM on top since they are thinest
lines(km,lwd=2,lty=c(1,2))



legend("topright",legend=c("6-MP","placebo"),lwd=2,lty=c(2,1))

#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenLeuk2Weibull.pdf")

s<-summary(sr)
confint(sr)

exp(coef(sr))
exp(confint(sr))


Beta<- summary(sr)$table[,"Value"]
Sigma<- sr$scale

HR<-exp(-Beta[2]/Sigma)
HR


vcov(sr)





## take derivative and do delta method
HRfunc<- function(Beta){
  exp( - Beta[2]*exp(-Beta[3]))
}

HRfunc(Beta)

HRci<-function(beta=Beta,srout=sr,conf.level=0.95,logtrans=TRUE){
     alpha<- 1-conf.level
     V<- vcov(srout)[2:3,2:3]
     HR<- HRfunc(beta)
     
    if (!logtrans){
      # Derivative function of Hazard Ratio
      dHR<- function(Beta){
        # get derivative=[dHR(b2,b3)/db2, dHR(b2,b3)/db3 ]
        out<- unname(HRfunc(Beta)* c( - exp(-Beta[3]), Beta[2]*exp(-Beta[3]) ) )
        names(out)<-c("dHR/dBetaTrt","dHR/dlogScale")
        out
      }
      # delta method
      D<- matrix(dHR(beta),1,2)
      VarHR<- as.numeric( D %*% V %*% t(D) )  
      ci<- unname(HR + c(-1,1)*qnorm(1-alpha/2)*sqrt(VarHR) )
    } else {
      dlogHR<-function(Beta){
        # get derivative of log(HR)
        out<- unname( c( - exp(-Beta[3]), Beta[2]*exp(-Beta[3]) ) )
        names(out)<-c("dlogHR/dBetaTrt","dlogHR/dlogScale")
        out
      }
      D<- matrix(dlogHR(beta),1,2)
      VarlogHR<- as.numeric( D %*% V %*% t(D) )  
      ci<- unname(exp(log(HR) + c(-1,1)*qnorm(1-alpha/2)*sqrt(VarlogHR) ))
    }
    out<-c(HR=HR,lower=ci[1],upper=ci[2])
    out
}


HRci(logtrans=FALSE)
HRci(logtrans=TRUE)

