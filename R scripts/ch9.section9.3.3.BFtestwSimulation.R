# Simulation at end of Section 9.3.3

source("ch9.bfTest.R")
source("ch9.BFdistn.R")

compareTConst<-function(r,n1,n2){
    dfWelch<- (r^2/n1 + 1/n2 )^2/( (r^4/(n1^2*(n1-1))) + (1/(n2^2*(n2-1)))  )
    c(n1=n1,n2=n2,w=qt(.975,dfWelch),bf=qbf(.975,n1,n2,s1=r,s2=1),h=qt(.975,min(c(n1,n2))),
      t=qt(.975,n1+n2-2))
}

nsim<-10^2
out<-matrix(NA,nsim,6,dimnames=list(NULL,c("n1","n2","w","bf","h","t")))
R<-.5
for (i in 1:nsim){
    N1<-rpois(1,5)+2
    N2<-rpois(1,5)+2
    s1sq<- rchisq(1,N1-1)
    s2sq<- R*rchisq(1,N2-1)    
    out[i,]<-compareTConst(sqrt(s1sq)/sqrt(s2sq),N1,N2)
}
out

pairs(cbind(pmin(out[,1],out[,2]),out[,-6]))
 plot(pmin(out[,1],out[,2]),out[,"h"] - out[,"bf"])




## simulate the distribution of the Behrens-Fisher Statistic
## Assume 
## Y1 ~ F1 = N(u1, v1)
## Y2 ~ F2 = N(u2,v2)
## take n1 from F1 and n2 from F2
## let ybar1 and ybar2 be means 
##   s1^2 and s2^2 be sample variances
## We want the distribution of:
##
##      (y2-y1) - (u2-u1)    
##  T= ------------------------
##      sqrt{s1^2/n1 +s2^2/n2}
##
##   Let rv = v1/v2

rbfs<-function(n1,n2,rv,nsim=10^4){
     Z2<-rnorm(nsim)
     Z1<- sqrt(rv)* rnorm(nsim)
     W2<-  (1/(n2*(n2-1)))*rchisq(nsim,n2-1)
     W1<-  (rv/(n1*(n1-1)))*rchisq(nsim,n1-1)
     T<- (Z2-Z1)/sqrt(W2+W1)
     T
}


hist(rbfs(10,4,2000))
 


Delta<-c(0,0,0,1,1,1)
Gamma<-c(.5,1,2,.5,1,2)
n2<-10
n1<- 20
nsim<-10^4

power<-matrix(NA,length(Delta),4)
propReject<-function(x){ length(x[x<=0.05])/length(x) }

alltests<-function(x,y){
  tout<-t.test(x,y, var.equal = FALSE)
  nmin<- min(length(x),length(y))
  pw<- tout$p.value
  ph<- 2*pt(-abs(tout$statistic),nmin-1)
  pt<- t.test(x,y, var.equal=TRUE)$p.value
  pbf<-bfTest(x,y)$p.value
  c(pt,pw,ph,pbf)
}

#alltests(1:5,3:12)
createData<-function(n,d,g){
  x<-rnorm(n[1])
  y<-rnorm(n[2],d,sd=sqrt(g))
  list(x=x,y=y)
}
#x<-createData(c(10^4,10^4),5,4)

set.seed(10301)
for (i in 1:length(Delta)){
  pvals<-matrix(NA,nsim,4)
  for (j in 1:nsim){
    x<-createData(c(n1,n2),Delta[i],Gamma[i])
    pvals[j,]<-alltests(x$x,x$y)
  }
  power[i,]<-apply(pvals,2,propReject)
}



powern1.20.n2.10<- cbind(Delta,Gamma,power)

dimnames(powern1.20.n2.10)<-list(NULL,c("Delta","v2/v1","t","Welch","Hsu","B-F"))

#write.csv(powern1.20.n2.10,file="ch2Ord.BFtestwSimulationResults.csv")



