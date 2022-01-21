source("ch9.bfTest.R")
source("ch9.BFdistn.R")


createData<-function(n,u=c(10,5),shift=c(0,5)){
  y1<-rpois(n[1],u[1])+ shift[1]
  y2<-rpois(n[2],u[2])+ shift[2]
  list(y1=y1,y2=y2)
}


nsim<-10^4
n<-c(20,10)
u<-c(10,5)
## Null hypothesis
#s<-c(0,5)
## Alternative hypothesis
s<-c(0,2.5)

pt<-pw<-pBF<-rep(NA,nsim)

set.seed(492191)
for (i in 1:nsim){
  x<-createData(n,u,s)
  pt[i]<-t.test(x$y1,x$y2, var.equal = TRUE)$p.value
  pw[i]<-t.test(x$y1,x$y2, var.equal = FALSE)$p.value
  pBF[i]<-bfTest(x$y1,x$y2)$p.value
}


propReject<-function(x){ length(x[x<=0.05])/length(x) }

propReject(pt)
propReject(pw)
propReject(pBF)
