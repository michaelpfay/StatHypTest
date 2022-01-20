# Fig 5.3
## plot ttestSimulation results (Power)
load("ch5.ttestSimulation1.RData")

par(mfrow=c(1,1),mar=c(6,4,4,2))


plot(c(0.5,4.5),c(0,100),type="n",xlab="",ylab="Percent Reject",
    main="Simulated Power",axes=FALSE)
axis(2)
bit<-.25
AT<-sort(c((1:4)-bit,1:4,(1:4)+bit))
LABLES<-rep(c(20,50,200),4)
axis(1,at=AT,label=LABLES,las=2)

Dnames<-c("Normal","Cauchy","t (df=2)","LogGamma(Shape=0.2)")
for (i in 1:3){
    mtext(Dnames[i],1,4,at=i)
}
mtext("Log-Gamma",1,3,at=4)
mtext("(Shape=.2)",1,4,at=4)

lines(c(-1,6),c(2.5,2.5),lwd=2,col="gray",lty=2)


## plot Power
PCH<-c(21:24)
#PCH<-c(15:18)
COL<-c("gray","black")

for (i in 1:4){

rejectLow<-c(
 sim.norm.20$rejectLow[i+1],
 sim.norm.50$rejectLow[i+1],
 sim.norm.200$rejectLow[i+1],
 sim.tdf1.20$rejectLow[i+1],
 sim.tdf1.50$rejectLow[i+1],
 sim.tdf1.200$rejectLow[i+1],
 sim.tdf2.20$rejectLow[i+1],
 sim.tdf2.50$rejectLow[i+1],
 sim.tdf2.200$rejectLow[i+1],
 sim.logGam.20$rejectLow[i+1],
 sim.logGam.50$rejectLow[i+1],
 sim.logGam.200$rejectLow[i+1])

rejectHi<-c(
 sim.norm.20$rejectHi[i+1],
 sim.norm.50$rejectHi[i+1],
 sim.norm.200$rejectHi[i+1],
 sim.tdf1.20$rejectHi[i+1],
 sim.tdf1.50$rejectHi[i+1],
 sim.tdf1.200$rejectHi[i+1],
 sim.tdf2.20$rejectHi[i+1],
 sim.tdf2.50$rejectHi[i+1],
 sim.tdf2.200$rejectHi[i+1],
 sim.logGam.20$rejectHi[i+1],
 sim.logGam.50$rejectHi[i+1],
 sim.logGam.200$rejectHi[i+1])


points(AT-bit/6,rejectLow,pch=PCH[i],col=COL[1])
points(AT+bit/6,rejectHi,pch=PCH[i],col=COL[2])
}

box()

#dev.print(pdf,"../graph/ch1OrdTtestSimP.pdf")