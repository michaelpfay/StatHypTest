# Fig 18.1
library(binseqtest)

#B<-designOBFpower(theta0=0.5,theta1=0.7,power=0.80, tsalpha=c(.025,.025))

B<-designOBF(30,k=3,theta0=0.5)
plot(B, rcol=c(black="black", gray="gray",black="black"))

#dev.print(pdf,file="../book/graph/chGS.OFk3.pdf")


sB<-stopTable(B)

dim(sB$table)


tab<-cbind(sB$table[,1:2],mle=sB$table[,"S"]/sB$table[,"N"],pL=sB$table[,"pL"])
tab<-tab[order(tab[,"mle"]),]
tab

powerBsb(B,theta=0.5)

prStop(B,theta=0.5)