# Figures 6.1 and 6.2
# takes a while to run (a couple of minutes)

seriesWSR<-function(x,y=NULL,theta=NULL,nt=100,package="stats",zmethod="Pratt"){
    if (is.null(y)){ 
        # if y is null, then x is the difference y-x
        d<-x
        # create fake x and y for input into certain packages
        x<-rep(0,length(d))
        y<-d
    } else { d<- y-x  }
    if (is.null(theta)){
        theta<- min(d) + ((1:nt)/(nt+1))*(max(d) - min(d) )
    }
    nt<-length(theta)
    pL<-pG<-rep(NA,nt)

    if (package=="stats"){
        ci<-wilcox.test(d,exact=TRUE,conf.int=TRUE)$conf.int
        for (i in 1:nt){
            dtemp<-d - theta[i]
            pL[i]<-wilcox.test(dtemp,exact=TRUE,conf.int=FALSE,alternative="less")$p.value
            pG[i]<-wilcox.test(dtemp,exact=TRUE,conf.int=FALSE,alternative="greater")$p.value
        }
    } else if (package=="exactRankTests"){
        ci<-wilcox.exact(d,exact=TRUE,conf.int=TRUE)$conf.int
        for (i in 1:nt){
            ytemp<-y -theta[i]
            pL[i]<-wilcox.exact(ytemp,x,paired=TRUE,conf.int=FALSE,alternative="less")$p.value
            pG[i]<-wilcox.exact(ytemp,x,paired=TRUE,conf.int=FALSE,alternative="greater")$p.value
        }
    } else if (package=="coin"){
        ci<-NULL
        for (i in 1:nt){
            ytemp<-y -theta[i]
            pL[i]<-pvalue(wilcoxsign_test(ytemp~x,zero.method=zmethod,alternative="less",distribution=exact()))
            pG[i]<-pvalue(wilcoxsign_test(ytemp~x,zero.method=zmethod,alternative="greater",distribution=exact()))
        }
    }
    list(theta=theta,pL=pL,pG=pG,ci=ci)
}



plotWSRseries<-function(sout,ci=NULL,YLIM=c(0,.1),plots="one.sided"){
    if (plots=="one.sided" | plots=="both"){
    plot(sout$theta,sout$pL, type="l",lwd=2,
        xlab=expression(theta[0]),col="gray",ylab="One-sided p-values", ylim=YLIM)
        lines(sout$theta,sout$pG, lwd=2)
    lines(range(sout$theta),c(.025,.025),lwd=1,lty=2)
    if (!is.null(sout$ci) | !is.null(ci)){
        if (is.null(ci)) ci<-sout$ci
        lines(c(ci[1],ci[1]),c(0,1))
        lines(c(ci[2],ci[2]),c(0,1))
    }
    }
    if (plots=="two.sided" | plots=="both"){
    pc<-pmin(1,2*sout$pL,2*sout$pG)
    plot(sout$theta,pc, type="l",lwd=2, xlab=expression(beta[0]),col="black",ylab="Two-sided p-values", ylim=c(0,1))
    lines(range(sout$theta),c(.05,.05),lwd=1,lty=2)
    }
}


d<-c(0,0,-2,-3,-5,6,9,11,12,15,16)
library(coin)

soutW<-seriesWSR(d,nt=2099,package="coin",zmethod="Wilcoxon")
soutP<-seriesWSR(d,nt=2099,package="coin",zmethod="Pratt")

# Fig 6.2
par(mfrow=c(2,1))
plotWSRseries(soutW,YLIM=c(0,.08))
title("Wilcoxon Method: Do Not Rank Zeros")
symbols(c(0,9),c(.029,.07),circles=c(.5,.5),inches=FALSE,col="gray",add=TRUE)
plotWSRseries(soutP,YLIM=c(0,.08))
title("Pratt Method: Rank Zeros")
#dev.print(pdf,file="../graph/chPairedWSRpvalues.pdf")

# Fig 6.1
plotWSRseries(soutW,YLIM=c(0,1))
title("Wilcoxon Method: Do Not Rank Zeros")
## Note all of the non-monotonic values, something at -2 I think
symbols(c(0,0,9,9),c(.97,.029,.07,.94),circles=c(.5,.5,.5,.5),inches=FALSE,col="gray",add=TRUE)

plotWSRseries(soutP,YLIM=c(0,1))
title("Pratt Method: Rank Zeros")
#dev.print(pdf,file="../graph/chPairedWSRpvalues2.pdf")

