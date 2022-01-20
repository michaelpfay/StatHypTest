# Section 7.4.4
# Calculations using Wang (2010) and Wang and Shan (2015) methods
install.packages("ExactCIdiff")
library(ExactCIdiff)


CL<-c(.95,.96)
CL2<-c(.975,.98)
CI<-CI2<-CI3<-CI4<-matrix(NA,length(CL),2)
for (i in 1:length(CL)){
  CI[i,]<-BinomCI(7,5,2,2,conf.level=CL[i],CItype="Two.sided")$ExactCI
  CI2[i,]<-BinomCI(5,7,2,2,conf.level=CL[i],CItype="Two.sided")$ExactCI
  CI3[i,]<-BinomCI(7,5,2,2,conf.level=CL2[i],CItype="Upper")$ExactCI
  CI4[i,]<-BinomCI(7,5,2,2,conf.level=CL2[i],CItype="Upper")$ExactCI
}

CI
CI2
CI3
CI4
BinomCI(7,5,2,2,conf.level=.95, CItype="Lower")

BinomCI(5,7,2,2,conf.level=.95, CItype="Lower")



BinomCI(10,10,1,9,conf.level=CL[i],CItype="Two.sided")$ExactCI