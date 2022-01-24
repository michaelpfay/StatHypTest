# Section 14.2.1, Problem 2
createData<-function(n){
  data.frame(y=rnorm(n),
             x1=rbinom(n,1,.5),
             x2=rbinom(n,1,.5),
             x3=rbinom(n,1,.5),
             x4=rbinom(n,1,.5),
             x5=rbinom(n,1,.5),
             x6=rnorm(n),
             x7=rnorm(n),
             x8=rnorm(n),
             x9=rnorm(n),
             x10=rnorm(n)
             )
}

set.seed(1)

nsim<-10^4
pvals<-matrix(NA,nsim,10)
for (i in 1:nsim){
  d<-createData(100)
  lout<-lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=d)
  pvals[i,]<-summary(lout)$coef[-1,4]
}

countReject<-function(p){ length(p[p<=0.05])}
nReject<- apply(pvals,1,countReject)
nReject2<- apply(pvals,2,countReject)
nReject2
table(nReject)