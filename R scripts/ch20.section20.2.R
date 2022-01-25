# Section 20.2
library(exactci)

N<- 10^(1:6)
CIL<- rep(NA,length(N))
for (i in 1:length(N)){
  CIL[i]<-powerBinom(n=N[i],type="cilength")$cilength
}

names(CIL)<-paste0("N=",N)
CIL

# It took a long time to calculate (a few minutes), so here are the results
#> CIL
#     N=10       N=100      N=1000     N=10000     N=1e+05     N=1e+06 
#0.599878392 0.202394744 0.062870719 0.019696259 0.006207844 0.001960961 