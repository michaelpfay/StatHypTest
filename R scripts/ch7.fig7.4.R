# Figure 7.4
# may take a while to run
library(exact2x2)

n1<-15
n2<-15

ciMatrix<-function(n1,n2,cifunc,...){
  x1<-0:n1
  x2<-0:n2
  lower<-upper<-matrix(NA,n1+1,n2+1,
                       dimnames=list(paste0("x1=",x1),paste0("x2=",x2)))
  for (i in 1:(n1+1)){
    for (j in 1:(n2+1)){
      ci<- cifunc(x1[i],n1,x2[j],n2,...)$conf.int
      lower[i,j]<-ci[1]
      upper[i,j]<-ci[2]
    }
  }
  list(lower=lower,upper=upper)
}

diffpropci<-function(x1,n1,x2,n2,conf.level=0.95){
  t2<-x2/n2
  t1<-x1/n1
  sigma<- sqrt( t1*(1-t1)/n1 + t2*(1-t2)/n2  )
  alpha<- 1- conf.level
  list(conf.int=t2-t1 + c(-1,1)*qnorm(1-alpha/2)*sigma,
       p.value= 2*(1-pnorm(abs((t2-t1)/sigma))))
}



cimUC<-ciMatrix(n1,n2,uncondExact2x2,conf.int=TRUE)
cimM<-ciMatrix(n1,n2,binomMeld.test)
cimMmp<-ciMatrix(n1,n2,binomMeld.test,midp=TRUE)
cimW<-ciMatrix(n1,n2,diffpropci)


coverage<-function(cim,thetas){
  nt<-nrow(thetas)
  d<-dim(cim$lower)
  n1<-d[1]-1
  n2<-d[2]-1
  errlo<-errhi<-cover<-rep(NA,nt)

  for (i in 1:nt){
    beta<- thetas[i,2] - thetas[i,1]
    f<- matrix(dbinom(0:n1,n1,thetas[i,1]),n1+1,1) %*% 
        matrix(dbinom(0:n2,n2,thetas[i,2]),1,n2+1)
    errlo[i]<- sum(f[beta<cim$lower])
    errhi[i]<- sum(f[beta>cim$upper])
    cover[i]<- sum(f[beta>=cim$lower & beta<=cim$upper])
  }
  list(thetas=thetas,errlo=errlo,errhi=errhi,cover=cover)
}



plotCoverage<-function(cov,title=""){
  plot(c(0,1.4),c(0,1),type="n", xlab="",
       ylab="",axes=FALSE,main="")
  axis(1,at=c(0:5/5),labels=c(0:5/5))
  axis(2,las=2)
  mtext(3,line=1,text=title,at=.5)
  mtext(1,line=3,text=expression(theta[1]),at=.5)
  mtext(2,line=3,text=expression(theta[2]),at=.5, las=2)
  #box()
  
  square<-function(t1,t2,b=1/200,COL){
    polygon(c(t1-b,t1-b,t1+b,t1+b,t1-b),
            c(t2-b,t2+b,t2+b,t2-b,t2-b),col=COL,border=NA)
  }
  
  nc<-dim(cov$thetas)[1]
  ## put error in integer units of 0.005
  ## 0 - 0.01 = 1
  ## 0.01- 0.02 = 2, etc
  cgrp<- cov$errlo/0.01
  ## put exact 0 a little higher so does not get numbered 0
  cgrp[cgrp<=0]<- 0.001
  cgrp<-ceiling(cgrp)
  cgrp[cgrp>6]<-6
  ## colors (gray(0)=black, gray(1)=white:

  clabel<-c("99-100","98-99","97-98","96-97","95-96","<95")
  COLS<-gray(c(.8,.9,1,.4,.3,0))
  
  for (i in 1:nc){
    square(cov$thetas[i,1],cov$thetas[i,2],COL=COLS[cgrp[i]])
  }
  
  for (i in 1:length(clabel)){
    square(1.1,i*.1,b=.04,COL=COLS[i])
    text(1.3,i*.1, labels=clabel[i])
  }
  
}



th<- 1:100/100 - 1/200
th1<- rep(th,length(th))
th2<-rep(th,each=length(th))
thetas<-cbind(th1,th2)

covUC<-coverage(cimUC,thetas)
covM<-coverage(cimM,thetas)
covMmp<-coverage(cimMmp,thetas)
covW<-coverage(cimW,thetas)

par(mfrow=c(2,2))

plotCoverage(covW,"Wald")
plotCoverage(covM,"Melded")
plotCoverage(covMmp,"mid-p Melded")
plotCoverage(covUC,"modified Boschloo")


#dev.print(pdf, file="../book/graph/ch2binCoveragePlots.pdf")
