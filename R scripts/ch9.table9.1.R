# Table 9.1
# Proportional odds to Mann-Whitney functional
## continuous case 

potomw<-function(psi,gamma=NULL){
  if (is.null(gamma)){
    # continuous case 
    out<- psi*(-1+psi-log(psi))/ (-1+psi)^2
    out[psi==1]<- 0.5
    out[psi==0]<-0
    out[psi==Inf]<- 1
    } else {
    # discrete case 
    n<-length(psi)
    out<-rep(NA,n)
    p<-length(gamma)
    for (i in 1:n){
      F1<- c(plogis(gamma),1)
      F2<- c(plogis(gamma-log(psi[i])),1)
      f1<- diff(c(0,F1))
      f2<- diff(c(0,F2))
      out[i]<- sum(F1[-(p+1)]*f2[-1]) + 0.5*sum(f1*f2)
    }
   }
  out
}


potomw(0)
potomw(0,c(-800:800)/100)

psi<-c(.25,.5,1,2,4)
phiL<-round(potomw(psi),3)
phi1<-round(potomw(psi,c(-1.1,-0.7,0.3,2.3)),3)
phi2<-round(potomw(psi,c(-1.5,.9)),3)



tab<-matrix(c(psi,phiL,phi1,phi2),nrow=length(psi),dimnames=list(NULL,c("psi","phiL","phi(1)","phi(2)")))

library(xtable)
print(xtable(tab,digits=c(0,2,3,3,3)),include.rownames=FALSE)


