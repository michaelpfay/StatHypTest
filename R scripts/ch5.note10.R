# Chapter 5 Note 10
## Simulation comparing simple mean of rates to quasi-Poisson

NSIM<-10^4
sim<-function(n=20,lambda=3,a=2,rtime=function(n){ 1+2*runif(n) },nsim=NSIM){
    m1<-m2<-m3<-m4<-rep(NA,nsim) 
    cilen1<-cilen2<-cilen3<-cilen4<-rep(NA,nsim)
    cov1<-cov2<-cov3<-cov4<-rep(1,nsim)
    for (i in 1:nsim){
        time<- rtime(n)
        #u<- time*lambda
        #a<- a/u
        w<-rgamma(n,1/a,scale=a)        
        x<-rpois(n,time*w*lambda)
        #x<-rpois(n,time*lambda)
        ci1<- t.test(x/time)$conf.int
        if (var(x/time)==0) ci1<-c(0,0)
        cilen1[i]<-ci1[2]-ci1[1]
        if (ci1[1]>lambda | ci1[2]<lambda) cov1[i]<-0
        m1[i]<-mean(x/time)
        g<-glm(x~1+offset(log(time)),family="quasipoisson")
        results<-summary(g)$coef
        beta<-results[1]
        seb<- results[2]
        ci2<- exp( beta+c(-1,1)*seb*qt(.975,n-1) )
        #ci2<-exp(confint(g))
        #ci2<-poisson.test(sum(x),sum(time))$conf.int
        cilen2[i]<-ci2[2]-ci2[1]
        if (ci2[1]>lambda | ci2[2]<lambda) cov2[i]<-0
        m2[i]<-exp( coef(g))
        # Fay-Feuer
        ff<-wspoisson.test(x,1/(n*time))
        m3[i]<-ff$estimate
        ci3<-ff$conf.int
        cilen3[i]<-ci3[2]-ci3[1]
        if (ci3[1]>lambda | ci3[2]<lambda) cov3[i]<-0
        #ff<-wspoisson.test(x,1/(n*time),midp=TRUE)
        #m4[i]<-ff$estimate
        #ci4<-ff$conf.int
        out4<-wspoisson.ci(x,time)
        m4[i]<-out4$estimate
        ci4<-c(out4$lower,out4$upper)
        cilen4[i]<-ci4[2]-ci4[1]
        if (ci4[1]>lambda | ci4[2]<lambda) cov4[i]<-0
   }

    c(n=n,lambda=lambda,a=a,MSE1=mean( (m1-lambda)^2 ),
         MSE2=mean( (m2-lambda)^2 ),
         MSE3=mean( (m3-lambda)^2 ),
         MSE4=mean( (m4-lambda)^2 ),
         cov1=mean(cov1),
         cov2=mean(cov2),
         cov3=mean(cov3),
         cov4=mean(cov4),
         cilen1=mean(cilen1),
         cilen2=mean(cilen2),
         cilen3=mean(cilen3),
         cilen4=mean(cilen4))
}


set.seed(29291)
N<-c(20,50,20,50,20,50)
A<-c(2,2,.8,.8,.2,.2)
LAMBDA<-c(rep(3,6),rep(10,6))
N<-rep(N,2)
A<-rep(A,2)

OUT<-matrix(NA,length(N),15)

for (i in 1:length(N)){
     outi<-sim(n=N[i],lambda=LAMBDA[i],a=A[i])
     OUT[i,]<-outi
}
dimnames(OUT)[[2]]<-names(outi)


#save(OUT,file="ch5.note10results.Rdata")
