library(exactci)

n<-10
x<-3
p<-.21
P<-c(1:999/1000)
P<-.21
OUT<-matrix(NA,length(P),3,dimnames=list(P,c("Blaker","MinLike","Central")))
for (i in 1:length(P)){
    OUT[i,]<-c(binom.exact(x,n,p=P[i],tsmethod="blaker")$p.value,
               binom.exact(x,n,p=P[i],tsmethod="minlike")$p.value,
               binom.exact(x,n,p=P[i],tsmethod="central")$p.value)
}
OUT
range(OUT[,1]-OUT[,2])

P[OUT[,1]-OUT[,2]< -.2]
P[OUT[,1]-OUT[,2]> .2]

# Figure 4.4
plotbinom<-function(x,n,p){
    d<-dbinom(0:n,n,p)
    plot(c(-1,n+1),range(c(0,d)),type="n",xlab="X",ylab="",bty="l")
    for (i in 1:(n+1)){
        b<-.5
        if (i==4){ COL<-"black"
        } else COL<-NA
        polygon(c(i-1-b,i-1+b,i-1+b,i-1-b,i-1-b),
            c(0,0,d[i],d[i],0),col=COL,border="black")
    }

}

plotbinom(x,n,.21)
dev.print(pdf,file="../graph/ch1bx3n10.pdf")

d<-dbinom(0:10,10,.21)
sum(d[4:11])
d[1]

#plotbinom(x,n,.36)
#binom.exact(x,n,p=.21,tsmethod="blaker")
#binom.exact(x,n,p=.21,tsmethod="minlike")