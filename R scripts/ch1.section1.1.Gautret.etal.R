# Recalculate Table 2 statistics from Gautret using different tests
library(exact2x2)
table2Stats<-function(x1,n1,x2,n2){
    I<- length(x1)
    pchisq<-pchisqC<-pfishC<-pfishML<-rep(NA,I)
    for (i in 1:I){
       m<- matrix(c(x1[i],n1[i]-x1[i],x2[i],n2[i]-x2[i]),2,2)
       pfishML[i]<-exact2x2(m, tsmethod="minlike")$p.value
       pfishC[i]<-exact2x2(m, tsmethod="central")$p.value
       pchisq[i]<-prop.test(c(x1[i],x2[i]),c(n1[i],n2[i]),correct=FALSE)$p.value
       pchisqC[i]<-prop.test(c(x1[i],x2[i]),c(n1[i],n2[i]),correct=TRUE)$p.value
    }
    out<-rbind(pfishML,pfishC,pchisq,pchisqC)
    dimnames(out)[[2]]<-paste0("day ",3:6)
    out
}

out<-table2Stats(c(10,12,13,14),c(20,20,20,20),c(1,4,3,2),c(16,16,16,16))
round(out,4)

