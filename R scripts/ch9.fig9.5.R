# Fig 9.5

par(mfrow=c(2,1))

## T-test with ARE=1 is df=18.76 (see fay and Proschan, 2010)
plot(c(-3,3),c(0,.5),type="n",axes=FALSE,xlab="",ylab="pdf",main="Scaled t distribution (df=18.76) vs standard normal")
box()
axis(1)
axis(2)

x<- -400:400/100

lines(x,dnorm(x),lwd=6,col=gray(.6))
lines(x,dt(x,df=18.76),lwd=3)

## plot gumble distribution
library(VGAM)


plot.gumbel<-function(loc,scale=sqrt(6)/pi){
  mu<- loc - scale*digamma(1)
  x<- (-500:500)/100 + mu
  y<-x-mu
  d<-dgumbel(y+mu,loc,scale)
  v<- pi^2 * scale^2/6
  lines(x-mu,d)
  out<-c(mean=mu,var=v)
  out
}

#plot(c(-5,5),c(0,1),type="n")

#plot.gumbel(3)
#plot.gumbel(1)


check.pdf<-function(pdf,x){
  DIFF<-diff(x)
  n<-length(x)
  d<- .5*pdf(x[-1]) + .5*pdf(x[-n]) 
  value<- .5*x[-1] + .5*x[-n]
  mu<-sum(DIFF*d*value)
  plot(x,pdf(x))
  out<-c(check=sum(DIFF*d),
         mean=mu, var=sum(DIFF*d*(value-mu)^2))
  out
}



pdfg<-function(x,a=1,b=10){
  dgamma(x,a,b)
}

pdf.loggamma<-function(x,a=1,b=10){
  (b^a/gamma(a))*exp(a*x - b*exp(x) )
}

pdf.sloggamma<-function(x,a=100,b=10){
  mu<-digamma(a)
  sigma<-sqrt(trigamma(a))
  z<- sigma*x + mu
  sigma*(b^a/gamma(a))*exp(a*z - b*exp(z) )
}


lines.loggamma<-function(a,b,x,LTY=1,...){
  pdf<-pdf.sloggamma(x,a=a,b=b)
  DIFF<-diff(x)
  n<-length(x)
  d<- .5*pdf[-1] + .5*pdf[-n] 
  value<- .5*x[-1] + .5*x[-n]
  mu<-sum(DIFF*d*value)
  pdf.plus.mu<-pdf.sloggamma(x+mu,a=a,b=b)
  #lines(x-mu,pdf.plus.mu)
  lines(x,pdf.plus.mu,lty=LTY,...)
  out<-c(check=sum(DIFF*d),
         mean=mu, var=sum(DIFF*d*(value-mu)^2))
  out
}
R<-(-500:500)/100
plot(c(-3,3),c(0,.5),type="n",axes=FALSE,xlab="",ylab="pdf",main="loggamma (a=5.55) vs standard normal")
box()
axis(1)
axis(2)

x<- -400:400/100

lines(x,dnorm(x),lwd=6,col=gray(.6))

#A<-c(.2,.5,1,5.55,100)
#LTYn<-c(6,4,5,1,2)
LTYn<-1
A<-c(5.55)
for (i in 1:length(A)){
  lines.loggamma(A[i],1,R,LTY=LTYn[i],lwd=3)
}


#dev.print(pdf,file="../graph/ch2OrdARE2plots.pdf")

#are.wt.loggamma<-function(a){
#  12*trigamma(a) * exp(2*lgamma(2*a) - 4*lgamma(a) - 4*a*log(2) )
#}
#ARE<-formatC(are.wt.loggamma(A),digits=3,format="f")
#space<-c("  ","  ","     ",""," ")
#LEG<-paste("a=",A,space," ARE=",ARE)
#LEG
#legend(c(-4),c(.8),LEG,lty=LTYn)


#dev.print(postscript,"H:\\main\\eps\\mpdr\\plot.loggamma.ps")
#dev.print(pdf,"H:\\My Documents\\methods\\mpdr\\tex\\plotloggamma.pdf")
#dev.print(pdf,"H:\\My Documents\\talks\\JSM 2008\\plotloggamma.pdf")