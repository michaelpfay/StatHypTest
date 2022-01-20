# Figure 5.2
## plot ttestSimulation results
load("ch5.ttestSimulation1.RData")

par(mfrow=c(1,1),mar=c(6,4,4,2))

## plot Type I error Rate
PCH<-c("L","G")
COL<-c("black","black")

plot(c(0.5,4.5),c(0,7.1),type="n",xlab="",ylab="Percent Reject",main="Simulated Type I Error Rate
Nominal=2.5%",axes=FALSE)
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

rejectLow<-c(
 sim.norm.20$rejectLow[1],
 sim.norm.50$rejectLow[1],
 sim.norm.200$rejectLow[1],
 sim.tdf1.20$rejectLow[1],
 sim.tdf1.50$rejectLow[1],
 sim.tdf1.200$rejectLow[1],
 sim.tdf2.20$rejectLow[1],
 sim.tdf2.50$rejectLow[1],
 sim.tdf2.200$rejectLow[1],
 sim.logGam.20$rejectLow[1],
 sim.logGam.50$rejectLow[1],
 sim.logGam.200$rejectLow[1])


rejectHi<-c(
 sim.norm.20$rejectHi[1],
 sim.norm.50$rejectHi[1],
 sim.norm.200$rejectHi[1],
 sim.tdf1.20$rejectHi[1],
 sim.tdf1.50$rejectHi[1],
 sim.tdf1.200$rejectHi[1],
 sim.tdf2.20$rejectHi[1],
 sim.tdf2.50$rejectHi[1],
 sim.tdf2.200$rejectHi[1],
 sim.logGam.20$rejectHi[1],
 sim.logGam.50$rejectHi[1],
 sim.logGam.200$rejectHi[1])


points(AT-bit/6,rejectLow,pch=PCH[1],col=COL[1])
points(AT+bit/6,rejectHi,pch=PCH[2],col=COL[2])

box()

#dev.print(pdf,"../graph/ch1OrdTtestSimSize.pdf")