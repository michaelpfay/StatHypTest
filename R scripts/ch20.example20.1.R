# example 20.1
# simulation
simttest<-function(nsim,n,delta,seed,varEqual=FALSE,sd2=1){
  set.seed(seed)
  # core function....7 lines
  p<-rep(NA,nsim)
  for (i in 1:nsim){
    y1<-rnorm(n)
    y2<-rnorm(n,mean=delta,sd=sd2)
    p[i]<-t.test(y1,y2,alternative="less",var.equal = varEqual)$p.value
  }
  length(p[p<=0.025])/nsim
  #
}

power<-rep(NA,10)
t0<-proc.time()
power[1]<-simttest(1e4,100,.5,1,varEqual=TRUE)
proc.time()-t0

for (i in 2:10){
  power[i]<-simttest(1e4,100,.5,i,varEqual=TRUE)
}
range(power)


powerW<-rep(NA,10)

for (i in 1:10){
  powerW[i]<-simttest(1e4,100,.5,i,varEqual=FALSE,sd2=1)
}
range(powerW)


powerW<-rep(NA,10)

for (i in 1:10){
  powerW[i]<-simttest(1e4,100,.5,i,varEqual=FALSE,sd2=2)
}
range(powerW)

# Simulate unequal variances
simttest(1e5,100,0,392013,varEqual=FALSE,sd2=2)

## Schouten (1999, Stat in MEd)
schouten<-function(alpha,gamma,r,xi,delta,sigma){
  ((r+xi)* sigma^2 * (qnorm(1-alpha)+qnorm(1-gamma))^2)/ (r*delta^2) +
    ((xi^2+r^3)*qnorm(1-alpha))/(2*r*(xi+r)^2)
}

schouten(.025,.4,1,4,.5,1)


#power.t.test(n=100,delta=delta, sig.level=.025,alternative="one.sided")$power





power<- .80
DELTA<- c(1:100/100)
calcSS<-function(power,DELTA){
  

normSS<-function(delta=.5,alpha=0.025,gamma=1-power){
  n<- (qnorm(1-alpha)+qnorm(1-gamma))^2  / delta^2
  n
}

normSS()


zn<-tn<-rep(NA,length(DELTA))

for (i in 1:length(DELTA)){
  zn[i]<- normSS(delta=DELTA[i])
  
  
  tn[i]<- power.t.test(delta=DELTA[i],sig.level=0.025,
               power=power,type="one.sample",
               alternative="one.sided")$n
}


names(zn)<-names(tn)<-paste0("Delta=",DELTA)

list(zn=zn,tn=tn)
}

x90<- calcSS(.9,DELTA)
x80<-calcSS(.8, DELTA)

range(x90$zn)
range(x80$zn)

range(x90$tn - x90$zn)
range(x80$tn - x80$zn)

plot(DELTA,x90$zn,type="l",col="gray",lty=1,lwd=6)
lines(DELTA,x90$tn,col="black",lty=1,lwd=2)

lines(DELTA,x80$zn,col=gray(.5),lty=1,lwd=6)
lines(DELTA,x80$tn,col=gray(0),lty=1,lwd=2)





