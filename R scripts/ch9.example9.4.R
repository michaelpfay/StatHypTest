# Example 9.4
#install.packages("mvtnorm")
library(mvtnorm)


## check to make sure we know what the output statistic from wilcox.test is
MWfunc<-function(x,y){
    n.x<-length(x)
    n.y<-length(y)
    r<-factor(rank(c(x,y)))
    rx<-r[1:n.x]
    ry<-r[-c(1:n.x)]
    fx<-table(rx)/n.x
    fy<-table(ry)/n.y
    Fx<- cumsum(fx)
    phiXY<- sum( (1-Fx)*fy + 0.5*fx*fy )
    ## compare to statistic from wilcox.test
    phiWilcox<- (1/(n.x*n.y))*(wilcox.test(x,y)$statistic)
    list(phiXY=phiXY,phiWilcox=phiWilcox)
}

set.seed(1)
x<-rpois(23,6)
y<-rpois(231,4)
MWfunc(x,y)



MWfuncBiNorm<-function(u1,v1,u2,v2,rho){
    pnorm( -(u2-u1)/sqrt(v1+v2-2*rho*sqrt(v1*v2)) )
}

#MWfuncBiNormOLD(1.2,3.1,2.1,4.7)
#MWfuncBiNorm(1.2,3.1,2.1,4.7,0)


simMixBiNorm<-function(n,meang,vg1,vg2,rhog,meanh,vh1,vh2,rhoh, lambda=0.5){
    ## part 1: simulate values 
    ng<-rbinom(1,n,lambda)
    nh<-n-ng
    vg12<- sqrt(vg1)*sqrt(vg2)*rhog
    sigmag<-matrix(c(vg1,vg12,vg12,vg2),2,2)
    vh12<- sqrt(vh1)*sqrt(vh2)*rhoh
    sigmah<-matrix(c(vh1,vh12,vh12,vh2),2,2)
    Yg<-rmvnorm(ng, meang, sigmag)
    Yh<-rmvnorm(nh, meanh, sigmah)
    Y<- rbind(Yg,Yh)
    phiXYest<-MWfunc(Y[,1],Y[,2])
    phiXYsim<-MWfunc(Y[,1],Y[,2])
    # fuction to estimate causal effect
    phiCausal<-function(Y){
        S<-sign(Y[,1] - Y[,2])
        n<-nrow(Y)
        (length(S[S==1]) + 0.5*length(S[S==0]))/n
    }
    phiCsim<-phiCausal(Y)
    ## Part 2: calculate values 
    ## assume Y ~ lambda*G + (1-lambda)*H
    phiXYint<- lambda^2 * MWfuncBiNorm(meang[1],vg1,meang[2],vg2,0) +
               lambda*(1-lambda)* MWfuncBiNorm(meanh[1],vh1,meang[2],vg2,0) +
               lambda*(1-lambda)* MWfuncBiNorm(meang[1],vg1,meanh[2],vh2,0) +
               (1-lambda)^2 * MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,0) 

    phiCint<- lambda*MWfuncBiNorm(meang[1],vg1,meang[2],vg2,rhog)+
              (1-lambda)*MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,rhoh)
    c(phiXYest=phiXYest,phiXYsim=phiXYsim,phiXYint=phiXYint,phiCsim=phiCsim,phiCint=phiCint)
}

## check the formula by simulation
set.seed(1)
simMixBiNorm(10^5,c(0,1),4,4,0,
             c(1,0),1,.1,.9,lambda=.7)
simMixBiNorm(10^5,c(2,3),1,1,.9,
             c(5,0),1,.1,0,lambda=.65)

MixBiNorm<-function(meang,vg1,vg2,rhog,meanh,vh1,vh2,rhoh, lambda=0.5){
   phiXYint<- lambda^2 * MWfuncBiNorm(meang[1],vg1,meang[2],vg2,0) +
               lambda*(1-lambda)* MWfuncBiNorm(meanh[1],vh1,meang[2],vg2,0) +
               lambda*(1-lambda)* MWfuncBiNorm(meang[1],vg1,meanh[2],vh2,0) +
               (1-lambda)^2 * MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,0) 
    phiCint<- lambda*MWfuncBiNorm(meang[1],vg1,meang[2],vg2,rhog)+
              (1-lambda)*MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,rhoh)
    nogood<- ((phiXYint<0.5) & (phiCint<0.5)) | ((phiXYint>0.5) & (phiCint>0.5))
    out<-c(phiXY=phiXYint,phiC=phiCint,
           stat=ifelse(nogood,0,abs(phiXYint-.5) + abs(phiCint-.5)))
    #    gC=MWfuncBiNorm(meang[1],vg1,meang[2],vg2,rhog)-.5,
    #    hC=MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,rhoh)-.5,
    #    gg=MWfuncBiNorm(meang[1],vg1,meang[2],vg2,0)-.5,
    #    hh=MWfuncBiNorm(meanh[1],vh1,meanh[2],vh2,0)-.5,
    #    gh=MWfuncBiNorm(meang[1],vg1,meanh[2],vh2,0)-.5,
    #    hg=MWfuncBiNorm(meanh[1],vh1,meang[2],vg2,0)-.5)
    out
}




wmwTestCalc<-function(x,y){
    n.x<-length(x)
    r<-rank(c(x,y))
    tab<-table(r)
    W<-sum(r[1:n.x]) - n.x*(n.x+1)/2

    Wwt<-wilcox.test(x,y)$statistic
    c(W1=W,W2=Wwt,tab=tab)
}



dg<- dh<- c(1,5)
rg<-rh<- c(.1,.5,1,1.5,2)
rhog<-rhoh<-c(0,.9)
lambda<-c(.65)
D<-c(1.5,2,2.5)

#dg<- dh<- c(1)
#rg<-rh<- c(1)
#rhog<-rhoh<-c(.9)
#lambda<-c(.9)
#D<-c(4)

OUT<-DNAME<-array(NA,dim=c(length(dg),length(dh),
                    length(rg),length(rh),
                    length(rhog),length(rhoh),
                    length(lambda), length(D)),
           dimnames=list(paste0("dg=",dg),
                         paste0("dh=",dh),
                         paste0("rg=",rg),
                         paste0("rh=",rh),
                         paste0("rhog=",rhog),
                         paste0("rhoh=",rhoh),
                         paste0("lambda=",lambda),
                         paste0("D=",D)))

prod(dim(OUT))

for (a in 1:length(dg)){
  for (b in 1:length(dh)){
    for (cc in 1:length(rg)){
      for (d in 1:length(rh)){
        for (e in 1:length(rhog)){
          for (f in 1:length(rhoh)){
            for (g in 1:length(lambda)){
              for (h in 1:length(D)){
               OUT[a,b,cc,d,e,f,g,h]<-MixBiNorm(c(0,dg[a])+D[h],1,rg[cc],rhog[e],
                        c(dh[b],0),1,rh[d],rhoh[f],lambda=lambda[g])["stat"]
               DNAME[a,b,cc,d,e,f,g,h]<-paste(a,b,cc,d,e,f,g,h,sep=",")
              }
            } 
          }
        }
      }
    }
  }
}


mix<-function(a,b,cc,d,e,f,g,h){
  parms<-paste("dg=",dg[a]," dh=",dh[b]," rg=",rg[cc]," rh=", rh[d], " rhog=", rhog[e],
        " rhoh=", rhoh[f], " lambda=", lambda[g], " D=", D[h])
  results<-MixBiNorm(c(0,dg[a])+D[h],1,rg[cc],rhog[e],
          c(dh[b],0),1,rh[d],rhoh[f],lambda=lambda[g])
  list(parms, results)
}
DNAME[OUT==max(OUT)]
mix(1,2,3,1,2,1,1,2)

MixBiNorm(c(0,1)+2,1,1,.9,c(5,0),1,.1,0,0.65)



#MixBiNorm(c(0,.1)+10,1,1,.9,
#                 c(.1,0),1,1,.9,lambda=.95)
#MixBiNorm(c(0,1),1,1,.9,
#                 c(5,0),1,1,.9,lambda=.7)




