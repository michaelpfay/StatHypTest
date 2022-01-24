# example 17.3 and Fig 17.1
# Do analysis method from Proschan et al 2001


# define p-values by Wald test
Zbin<-function(p1,p2,n1,n2){
  p0<- (n1*p1 + n2*p2)/(n1+n2) 
  Z<- (p1-p2)/sqrt( p0*(1-p0)*(1/n1 + 1/n2) )
  list(Z=Z,p.value=2*(1-pnorm(abs(Z))))
}
#Zbin<-function(p1,p2,n1,n2){
#  prop.test(round(c(p1*n1,p2*n2)),n=c(n1,n2),correct=TRUE)$p.value
#}



# observed case
Zbin(75/193,105/192,193,192)
# worst case
Zbin((75+11)/204,105/212,204,212)
# opposite arm imputation
Zbin((75+11*(105/192))/204,(105+20*(75/193))/212,204,212)
# pooled imputation
pp<-  (75+105)/(193+192)
Zbin((75+11*pp)/204,(105+20*pp)/212,204,212)



library(exact2x2)
library(TippingPoint)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(bayesSurv)
# Now consider testing by 
# 

x1<-75
n1<-204
x2<-105
n2<-212
m1<-11
m2<-20
treat<- c(rep(1,n1),rep(0,n2))
y<- c(rep(1,x1),rep(0,n1-m1-x1),rep(NA,m1),
      rep(1,x2),rep(0,n2-m2-x2),rep(NA,m2)
      )

# modify the TippingPoint R package so that if calculates central Fisher's exact p-values 
# instead of prop.test default p-values for binary data
source("ch17.TippingPoint.Modified.R")
class(y)<-"Modified"

# Fig 17.1
TippingPoint.modified(outcome=y, treat= treat,
             group.infor=TRUE, plot.type = "p.value",ind.values = TRUE,
             impValuesT  = NA,  impValuesColor = ,
             summary.type = "density", alpha = 0.95)



packageVersion("TippingPoint")

#dev.print(pdf, file="../book/graph/chMiBinaryMNAR1.pdf")


# Below is a kind of ugly version 
calcBinaryPvalues<-function(x1,n1,x2,n2,m1,m2, type="centalFishersExact",...){
  getpval<-function(X1,N1,X2,N2){
    fisher.exact(matrix(c(X1,N1-X1,X2,N2-X2),2,2),tsmethod="central")$p.value
  }
  #out<-matrix(NA,m1+1,m2+1,dimnames=list(paste0("xm1=",0:m1),paste0("xm2=",0:m2)))
  out<-matrix(NA,m1+1,m2+1,dimnames=list(0:m1,0:m2))
  xm1<- 0:m1
  xm2<- 0:m2
  
  for (i in 1:length(xm1)){
    for (j in 1:length(xm2)){
      out[i,j]<-getpval(x1+xm1[i],n1,x2+xm2[j],n2)
      # to check output is correct
      #out[i,j]<- xm1[i]+xm2[j]/100
    }
  }
  out
  #plot(0:m1, 0:m2,...)
  
}

z<-calcBinaryPvalues(75,204,105,212,11,20)
m1<-11
m2<-20
d<- expand.grid(x=0:m1,y=0:m2)
dim(d)
d$p.value<-as.vector(z)
ggplot(d,aes(x,y,fill=p.value)) + geom_raster()
