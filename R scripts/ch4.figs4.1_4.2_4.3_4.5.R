n<-30
L<-U<-Lbl<-Ubl<-Lmidp<-Umidp<-Lwilson<-Uwilson<-rep(NA,n+1)
library(BlakerCI)
library(PropCIs)

library(exactci)
## not done with making package yet, just use functions 
#DIR<-"H:/My Documents/methods/PUBLISHED/2010 Fay R Journal/exactci/R/"
#source(paste0(DIR,"binom.exact.R"))
#source(paste0(DIR,"binomControl.R"))
#source(paste0(DIR,"exactbinomCI.R"))
#source(paste0(DIR,"exactbinomPlot.R"))
#source(paste0(DIR,"exactbinomPvals.R"))


conflev<-0.95
for (i in 1:(n+1)){
    ci<-binom.test(i-1,n,conf.level=conflev)$conf.int
    L[i]<-ci[1]
    U[i]<-ci[2]
    ci2<-binom.blaker.limits(i-1, n, level = conflev, tol = 1e-10)
    Lbl[i]<-ci2[1]
    Ubl[i]<-ci2[2]
    ci3<-binom.exact(i-1,n,conf.level=conflev,midp=TRUE)$conf.int
    Lmidp[i]<-ci3[1]
    Umidp[i]<-ci3[2]
    
    ci4<- scoreci(i-1,n,conf.level=conflev)$conf.int
    Lwilson[i]<-ci4[1]
    Uwilson[i]<-ci4[2]
    
}



P<- (0:10000)/10000
errlo<-errhi<-errloBL<-errhiBL<-errloMidp<-errhiMidp<-errloW<-errhiW<-rep(NA,length(P))

for (j in 1:length(P)){
    d<-dbinom(0:n,n,P[j])
    errlo[j]<-sum(d[P[j]<L])
    errhi[j]<-sum(d[U<P[j]])
    errloBL[j]<-sum(d[P[j]<Lbl])
    errhiBL[j]<-sum(d[Ubl<P[j]])

    errloMidp[j]<-sum(d[P[j]<Lmidp])
    errhiMidp[j]<-sum(d[Umidp<P[j]])

    errloW[j]<-sum(d[P[j]<Lwilson])
    errhiW[j]<-sum(d[Uwilson<P[j]])
}

## Figure 4.3

par(mfrow=c(1,1))
plot(rep(0:n,2),c(L,U),type="n",xlab="x",ylab="")
#segments(0:n,Lbl,0:n,Ubl,lwd=16,col="gray")
for (i in 0:n){
    bit<-.2
    polygon(c(i-bit,i+bit,i+bit,i-bit,i-bit),
        c(Lbl[i+1],Lbl[i+1],Ubl[i+1],Ubl[i+1],Lbl[i+1]),
        border=NA,col="gray")
    bit<-.05
    polygon(c(i-bit,i+bit,i+bit,i-bit,i-bit),
        c(L[i+1],L[i+1],U[i+1],U[i+1],L[i+1]),
        border=NA,col="black")

}

## percent of exact central
## 
round(100-100*(Ubl-Lbl)/(U-L),1)
mean(100-100*(Ubl-Lbl)/(U-L))



#segments(0:n,L,0:n,U,lwd=3)
dev.print(pdf,file="../graph/ch1bCI30.pdf")

# Figure 4.1
#par(mfrow=c(1,3))
par(mfrow=c(3,1),mar=c(1.5,4,1,0.5),oma=c(2,0,2,0)+.5)
plot(P,errlo+errhi,type="l",lwd=1,ylim=c(0,.06),ylab="Total Error",xaxs="r")
lines(c(-1,10),c(0.05,0.05),lty=2,col=gray(.8))
plot(P,errhi,type="l",lwd=1,ylim=c(0,.05),ylab="Upper Limit Error",xaxs="r")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
plot(P,errlo,type="l",lwd=1,ylim=c(0,.05),ylab="Lower Limit Error",xaxs="r")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
mtext(paste0("Error for Clopper-Pearson 95% Confidence Interval (n=",n,")"),side=3,outer=TRUE)
mtext(expression(theta),side=1,line=1.2,outer=TRUE)


dev.print(pdf,file="../graph/ch1bCIerr1.pdf")
#dev.print(postscript,file="../graph/ch1bCIerr1.ps")

# Figure 4.5
par(mfrow=c(3,1),mar=c(1.5,4,1,0.5),oma=c(2,0,2,0)+.5)
plot(P,errloMidp+errhiMidp,type="l",lwd=1,ylim=c(0,.08),ylab="Total Error")
lines(c(-1,10),c(0.05,0.05),lty=2,col=gray(.8))
plot(P,errhiMidp,type="l",lwd=1,ylim=c(0,.08),ylab="Upper Limit Error")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
plot(P,errloMidp,type="l",lwd=1,ylim=c(0,.08),ylab="Lower Limit Error")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
mtext(paste0("Error for Mid-p 95% Confidence Interval (n=",n,")"),side=3,outer=TRUE)


dev.print(pdf,file="../graph/ch1bCIerrMidp.pdf")


par(mfrow=c(3,1),mar=c(1.5,4,1,0.5),oma=c(2,0,2,0)+.5)
plot(P,errloW+errhiW,type="l",lwd=1,ylim=c(0,.08),ylab="Total Error")
lines(c(-1,10),c(0.05,0.05),lty=2,col=gray(.8))
plot(P,errhiW,type="l",lwd=1,ylim=c(0,.08),ylab="Upper Limit Error")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
plot(P,errloW,type="l",lwd=1,ylim=c(0,.08),ylab="Lower Limit Error")
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
mtext(paste0("Error for Wilson 95% Confidence Interval (n=",n,")"),side=3,outer=TRUE)


dev.print(pdf,file="../book/graph/ch1bCIerrWilson.pdf")




# Figure 4.2
par(mfrow=c(3,1))
plot(P,errloBL+errhiBL,type="l",lwd=1,ylim=c(0,.06),ylab="Total Error",xlab=expression(theta))
lines(c(-1,10),c(0.05,0.05),lty=2,col=gray(.8))
plot(P,errhiBL,type="l",lwd=1,ylim=c(0,.06),ylab="Upper Limit Error",xlab=expression(theta))
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
plot(P,errloBL,type="l",lwd=1,ylim=c(0,.06),ylab="Lower Limit Error",xlab=expression(theta))
lines(c(-1,10),c(0.025,0.025),lty=2,col=gray(.8))
mtext(paste0("Error for Blaker 95% Confidence Interval (n=",n,")"),side=3,outer=TRUE)
mtext(expression(theta),side=1,line=1.2,outer=TRUE)


dev.print(pdf,file="../graph/ch1bCIerrBL.pdf")


