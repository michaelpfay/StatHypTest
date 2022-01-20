# Ch 6 Figures 6.3 and 6.4



# I got this data from 
# H:/My Documents/malaria/phase1/ama1/ama1cpg Mali adults/elisa/R/AMA1CPGMaliELISA.txt.R"
# I ran that, and saw the code that 
# reproduced Fig 1 of Sagara, et al 2009 Vaccine 27:7292-98.
# it called another file, H:/.../plot.elisa.mean.noSubjLines.R
# I stuck a broswer into that file and saved the day 0 data and the day 210 data 
# for the no CPG group.
# Then I saved the data
load("ch6elisa.RData")

gmeanCI<-matrix(NA,2,3,dimnames=list(c("day 0","day 210"),
     c("Gmean","Lower 95CL","Upper 95CL")))
t0<-t.test(log10(DFig1[,"avgELISA0"]))
gmeanCI[1,1]<-10^t0$estimate
gmeanCI[1,2:3]<-10^t0$conf.int
t0<-t.test(log10(DFig1[,"avgELISA210"]))
gmeanCI[2,1]<-10^t0$estimate
gmeanCI[2,2:3]<-10^t0$conf.int
gmeanCI

# Fig 6.3
par(mfrow=c(1,2))
d0<-round(DFig1[,"avgELISA0"],2)
d210<-round(DFig1[,"avgELISA210"],2)
t.test(d0,d210)
t.test(log10(d0),log10(d210))
paste0("(",d0,",",d210,")")

plot(c(rep(0,12),rep(210,12)),c(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"]),xlab="Vaccine Day",
   ylab="AMA1-C1 ELISA (ug/mL)",xlim=c(-100,310),axes=FALSE)
segments(rep(0,12),DFig1[,"avgELISA0"],rep(210,12),DFig1[,"avgELISA210"])
axis(1,at=c(0,210),labels=c(0,210))
axis(2)
box()
title("(a)")

plot(c(rep(0,12),rep(210,12)),c(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"]),xlab="Vaccine Day",
   log="y",ylab="AMA1-C1 ELISA (ug/mL)",xlim=c(-100,310),axes=FALSE)
segments(rep(0,12),DFig1[,"avgELISA0"],rep(210,12),DFig1[,"avgELISA210"])
axis(1,at=c(0,210),labels=c(0,210))
axis(2)
box()
title("(b)")
#dev.print(pdf,"../graph/chPaiElisaApplication.pdf")


# Fig 6.4
par(mfrow=c(2,1))
z<-rep(0,12)
plot(d210-d0,z,ylab="",xlab="ELISA(Day 210) - ELISA(Day 0)",axes=FALSE)
axis(1)
box()
title("(a) Difference")
plot(log(d210)-log(d0),z,log="",ylab="",xlab="ELISA(Day 210)/ELISA(Day 0)",axes=FALSE)
AT<-log(c(1,2,4,8,16,32))
axis(1,at=AT,labels=round(exp(AT),1))
box()
title("(b) Ratio (log scale)")
#dev.print(pdf,"../graph/chPaiElisaApplication2.pdf")



## arithmetic scale
t.test(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"])
t.test(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"],paired=TRUE)
cor.test(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"])
plot(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"],log="xy")
# log scale
t1<-t.test(log(DFig1[,"avgELISA210"]),log(DFig1[,"avgELISA0"]))
t1
exp(t1$conf.int)
exp(t1$estimate)
t2<-t.test(log(DFig1[,"avgELISA210"]),log(DFig1[,"avgELISA0"]),paired=TRUE)
t2
exp(t2$conf.int)
exp(t2$estimate)
## percent change 
100*(exp(t2$estimate)-1)
100*(exp(t2$conf.int)-1)

library(asht)
w1<-wsrTest(log(DFig1[,"avgELISA210"]),log(DFig1[,"avgELISA0"]))
w1
exp(w1$estimate)
exp(w1$conf.int)


t.test(log(DFig1[,"avgELISA0"]),log(DFig1[,"avgELISA210"]),paired=TRUE)
cor.test(DFig1[,"avgELISA0"],DFig1[,"avgELISA210"])
