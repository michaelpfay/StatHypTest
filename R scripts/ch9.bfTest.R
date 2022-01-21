# Behrens-Fisher test
# now available in asht R package
bfTest<-function(x,y,nsim=10^5,beta=0,conf.level=0.95,calcMethod=c("int","mc")){
    ## Behrens-Fisher test
    u2<-mean(x)
    s2<-sqrt(var(x))
    n2<-length(x)
    u1<-mean(y)
    s1<-sqrt(var(y))
    n1<-length(y)
    if (n1<2 | n2<2) stop("n1 and n2 must be at least 2")
    alpha<-1-conf.level   

    calcMethod<-match.arg(calcMethod)
    
    cbf<-NA
    if (calcMethod=="mc"){
       R2<- (u2 + (s2/sqrt(n2))*rt(nsim,n2-1) )  
       R1<- (u1 + (s1/sqrt(n1))*rt(nsim,n1-1) )
       D<- R2-R1
       CI<-quantile(D,probs=c(alpha/2,1-alpha/2))
       p.L<- length(D[D<=beta])/nsim
    } else if (calcMethod=="int"){
        ifunc<- function(x,b=beta){
            pt( (b+u1-u2+(s1/sqrt(n1))*x)/(s2/sqrt(n2)),n2-1) * dt(x,n1-1)
        }
        p.L<-integrate(ifunc,-Inf,Inf)$value
        cbf<-qbf(1-(1-conf.level)/2,n1,n2,s1=s1,s2=s2)
        CI<- (u2-u1) + c(-1,1)*cbf*sqrt(s1^2/n1 + s2^2/n2 )
        #CI<-c(NA,NA)
    }
    p.U<-1-p.L
    p<- min(2*p.L,2*p.U)
    list(t=u2-u1,p.L=p.L,p.U=p.U,p.value=p,cbf=cbf,cilo=CI[1],cihi=CI[2],conf.int=CI)
}

#set.seed(1321)
#x<-rnorm(10)
#y<-rnorm(4)+.01

#bfTest(x,y)
#bfTest(x,y,nsim=10^5, calcMethod = "mc")
#t.test(x,y)$p.value