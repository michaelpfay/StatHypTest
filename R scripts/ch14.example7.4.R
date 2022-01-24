# Ch 14, Problem P6, Example 7.4 
y<-c(4,0)
m<-c(20,24)
g<-c("SOC","ZMapp+SOC")

gout<-glm(cbind(y,m-y)~g,family=binomial())
summary(gout)