# Fig 21.1 and example 21.2

x<-12
n<-63

par(mfrow=c(2,1))

p<-0:1000/1000

thetaPrior<-0.075
b<-20
a<- round(thetaPrior*b/(1-thetaPrior),2)
a/(a+b)
a

#a<-b<-.5
plot(p, dbeta(p, a,b),type="l",lwd=4,xlab=expression(theta),ylab="",main=paste0("Prior: Beta(",a,",",b,")") )


plot(p, dbeta(p, a+x,b+n-x),type="l",lwd=4,xlab=expression(theta),ylab="",
      main=paste0("Posterior: Beta(",round(a+x,2),",",round(b+n-x,2),")") )

lines(c(x/n,x/n),c(0,15),col=gray(.7),lwd=4)
lines(p,dbeta(p, a+x,b+n-x),lwd=4)

#dev.print(pdf,file="../graph/chBTPrev.pdf")

## Calculation of some Bayesian results:

meanPosterior<- (a+x)/(n+a+b)
meanPosterior
modePosterior<- (a+x-1)/(n+a+b-2)
modePosterior


### Credible intervals.....

## Prior based on a and b
qbeta(0.025, a+x,b+n-x)
qbeta(0.975, a+x,b+n-x)


## Jeffery's prior 
qbeta(0.025, .5+x,.5+n-x)
qbeta(0.975, .5+x,.5+n-x)

binom.test(x,n)$conf.int

### Large sample size
x<-1200
n<-6300

## Prior based on a and b
qbeta(0.025, a+x,b+n-x)
qbeta(0.975, a+x,b+n-x)



## Jeffery's prior 
qbeta(0.025, .5+x,.5+n-x)
qbeta(0.975, .5+x,.5+n-x)

binom.test(x,n)$conf.int

binom.test(52263471,104490000)

