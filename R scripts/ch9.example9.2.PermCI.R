# Example 9.2 permutation CI part
# calculate permutation p-value and shift CI using difference in means permutation test
library(perm)

y1<-c(997, 1026,  953,  967,  996 )
y2<-c(1112, 1115, 1113, 1090, 1073)

cm10c5<-chooseMatrix(10,5)

rootfunc<-function(delta,alpha=0.025,alt="less"){
  permTS(y1,y2-delta, method="exact.ce",alternative=alt,control=permControl(cm=cm10c5))$p.value - alpha
}

rootfunc(0,alt="greater")
rootfunc(10^3,alt="greater")


uniroot(rootfunc,c(0,10^3),alt="less")$root

uniroot(rootfunc,c(0,10^3),alt="greater")$root


## PCLT approx
v1<- mean( (y1 - mean(y1) )^2 )
v2<- mean( (y2 - mean(y2) )^2 )
n1<- length(y1)
n2<-length(y2)
n<-n1+n2
q<-qnorm(0.975)
D<- mean(y2) - mean(y1)
CId<- q*sqrt((n/(n-q^2-1))*(v1/n2 + v2/n1))
lower<- D - CId
upper<- D + CId
lower
upper 
