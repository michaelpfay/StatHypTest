# Figure 7.3
# This requires a long time to run, over 5 hours originally


library(exact2x2)

ni<-10^c(1:5)
x1<- round(.8*ni)
x2<-round(.9*ni)

CIWald<-CIMeld<-CIUncond<-CImidpMeld<-matrix(NA,length(ni),2)
pWald<-pMeld<-pUncond<-pmidpMeld<-rep(NA,length(ni))

diffpropci<-function(x1,n1,x2,n2){
  t2<-x2/n2
  t1<-x1/n1
  sigma<- sqrt( t1*(1-t1)/n1 + t2*(1-t2)/n2  )
  list(conf.int=t2-t1 + c(-1,1)*qnorm(.975)*sigma,
       p.value= 2*(1-pnorm(abs((t2-t1)/sigma))))
}


t0<-proc.time()

for (i in 1:length(ni)){
  s<-diffpropci(x1[i],ni[i],x2[i],ni[i])
  pWald[i]<-s$p.value
  CIWald[i,]<-s$conf.int
}

t1<-proc.time()

for (i in 1:length(ni)){
   s<-binomMeld.test(x1[i],ni[i],x2[i],ni[i],parmtype="difference")
   pMeld[i]<-s$p.value
   CIMeld[i,]<-s$conf.int
}

t2<-proc.time()


for (i in 1:length(ni)){
  s<-binomMeld.test(x1[i],ni[i],x2[i],ni[i],parmtype="difference",midp=TRUE)
  pmidpMeld[i]<-s$p.value
  CImidpMeld[i,]<-s$conf.int
}

t3<-proc.time()

## this takes a long time (5 hours or so?)....
for (i in 1:3){
  s<-uncondExact2x2(x1[i],ni[i],x2[i],ni[i],parmtype="difference", conf.int=TRUE)
  pUncond[i]<-s$p.value
  CIUncond[i,]<-s$conf.int
}

t4<-proc.time()


timeWald<- t1-t0
timeMeld<- t2-t1
timemidpMeld<- t3-t2
timeUncond<- t4-t3
time1<-c(timeWald,timeMeld,timeUncond)

CI1<-cbind(CIWald,CIMeld,CIUncond)

p1<-cbind(pWald,pMeld,pUncond)

time2<-c(timeWald,timeMeld,timemidpMeld)
CI2<-cbind(CIWald,CIMeld,CImidpMeld)
p2<-cbind(pWald,pMeld,pmidpMeld)

#save(time1,CI1,p1,file="ch2binWaldMeldUncond2.Rdata")
#load("ch2binWaldMeldUncond2.Rdata")

par(mfrow=c(1,1))

plot(c(min(CI2),max(CI2)),c(min(ni),max(ni)*3),log="y",xlab=expression(beta[d]),ylab=expression(n[a]),type="n")
COL<-c("gray","gray","black","black")
LTY<-c(1,3,1,3)
LWD<-c(8,5,3,1)
for (i in 1:3){
  lines(CI2[,2*i-1],ni,col=COL[i],lty=LTY[i],lwd=LWD[i])
  lines(CI2[,2*i],ni,col=COL[i],lty=LTY[i],lwd=LWD[i])
}

for (i in 3:3){
  lines(CI1[1:3,2*i-1],ni[1:3],col=COL[4],lty=LTY[4],lwd=LWD[4])
  lines(CI1[1:3,2*i],ni[1:3],col=COL[4],lty=LTY[4],lwd=LWD[4])
}

o<-c(1,3,4,2)
legend(.22,10^5,legend=c("Wald","Melded","mid-p Melded","modified Boschloo")[o],lty=LTY[o],lwd=LWD[o],col=COL[o])

#dev.print(pdf,file="../book/graph/ch2binWaldMeldUncond.pdf")
