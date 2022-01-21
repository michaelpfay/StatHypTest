# Section 12.2.4


## the idea is to compare the difference in means (ignoring cluster) 
## to difference in means of cluster means
## with random effects difference in means
##
## I want to show that the random effects picks out the right one,
## when the varw>>varb then ignore cluster, 
## and when varw<<varb use difference in mean of cluster means
## in between cases uses the random effects model. 

t0<-proc.time()

## Run same simulation as chStr.simRandomEffectsCluster.R (i.e., ch12.section12.2.3.R)
## Except: 1. make errors logistic
#          2. round responses to nearest integer
#          2. Use ordinal R package with and without random effects for cluster 
#library(lme4)
#library(lmerTest)
library(ordinal)
library(polr)
library(MASS)
nsim<- 10^4
set.seed(104301)

simCluster<-function(N1=10,N2=10,mu1=0,mu2=1, mClustSet=c(10,20,40), varw=1, varb=1){
   N<-N1+N2
   m<- sample(mClustSet,N, replace=TRUE)
   # b<- rnorm(N,sd=sqrt(varb))
   b<- rlogis(N,scale=sqrt(varb)*sqrt(3)/pi)
   B<- rep(b,m)
   CID<- factor(rep(1:N,m))
   MU<- rep(c(mu1,mu2),c(sum(m[1:N1]),sum(m[(N1+1):N])))
   G<- factor(rep(c("g1","g2"),c(sum(m[1:N1]),sum(m[(N1+1):N]))))
   #Y<- MU+B+ rnorm(sum(m),sd=sqrt(varw))
   Y<- MU+B+ rlogis(sum(m),scale=sqrt(varw)*sqrt(3)/pi)
   data.frame(MU,B,Y,G,CID,Ycat=factor(round(Y)))
}



d<-simCluster()


#x<-polr(Ycat~G,data=d, Hess=TRUE)


analyzeOrd<-function(d,trueDelta=0){
   # one obs one vote

   c1<- clmm2(Ycat~G, data=d, Hess=TRUE)
   s1<-summary(c1)$coefficients["Gg2",]
   ci<- unname(s1["Estimate"] + c(-1,1)*s1["Std. Error"]*qnorm(.975) )
   t1<-c(estimate=unname(s1["Estimate"]),
         lo=ci[1], hi=ci[2],
         p.value=unname(s1["Pr(>|z|)"]))

   c2<- clmm2(Ycat~G, random=CID, data=d, Hess=TRUE)
   s2<-summary(c2)$coefficients["Gg2",]
   ci<- unname(s2["Estimate"] + c(-1,1)*s2["Std. Error"]*qnorm(.975) )
   t2<-c(estimate=unname(s2["Estimate"]),
         lo=ci[1], hi=ci[2],
         p.value=unname(s2["Pr(>|z|)"]))
   
   
      extraStats<- function(s,trueDelta=1){
      c(s, cilen=unname(s["hi"]-s["lo"]),
        cover=unname(ifelse(s["lo"]<=trueDelta & s["hi"]>=trueDelta,1,0)),
        reject=unname(ifelse(s["p.value"]<=0.05,1,0)))
   }
   t1<-extraStats(t1,trueDelta)
   t2<-extraStats(t2,trueDelta)
   list(t1=t1,t2=t2)
}

temp<-analyzeOrd(d,1)
VARB<- c(.01,.1,1,10)^2
VARB<-c(.01,.1,1)^2
out<-array(NA,c(2,7,nsim,length(VARB)),dimnames=list(c("prop odds","prop odds RE"),
                                                       names(temp$t1),
                                                       1:nsim,
                                                     paste0("varb=",VARB)
                                                     ))
for (j in 1:length(VARB)){
   for (i in 1:nsim){
      d<-simCluster(mu2=0,varb=VARB[j])
      # trueDelta= mu2-mu1
      temp<- analyzeOrd(d,trueDelta = 0)
      out[1, ,i,j]<-temp$t1
      out[2, ,i,j]<-temp$t2
      print(paste("j=",j," i=",i))
   }
}
   

summarizeSim<-function(out,Delta=1,NA.RM=TRUE){

   getmean<-function(x){ mean(x["estimate",],na.rm=NA.RM) }
   getMSE<- function(x){ mean((x["estimate",]-Delta)^2,na.rm=NA.RM) }
   getCover<- function(x){ mean(x["cover",],na.rm=NA.RM)}
   getPower<- function(x){ mean(x["reject",],na.rm=NA.RM)}
   getCIlen<- function(x){ mean(x["cilen",],na.rm=NA.RM)}
   sout<- array(NA,c(2,length(VARB),5),
                dimnames=list(c("prop odds","prop odds RE"),
                              paste0("varb=",VARB),
                              c("mean","MSE","coverage","power","cilen")))
   sout[ , ,"mean"]<- apply(out,c(1,4),getmean)
   sout[ , ,"MSE"]<- apply(out,c(1,4),getMSE)
   sout[ , ,"coverage"]<- apply(out,c(1,4),getCover)
   sout[ , ,"power"]<- apply(out,c(1,4),getPower)
   sout[ , ,"cilen"]<- apply(out,c(1,4),getCIlen)
   sout
}

simSum<-summarizeSim(out,Delta=0)

#save(simSum, file="chStrsimRandomEffectsClusterOrdinal.RData")
#load("chStrsimRandomEffectsClusterOrdinal.RData")


plotSimSum<-function(ss,TITLE="",legloc=NULL){
   rclust<- ss[2,]
   rind<- ss[1,]
   plot(range(sqrt(VARB)),range(c(rclust,rind)),type="n",log="xy",xlab=expression(sigma[b]),ylab="ratio over lme",main=TITLE)
    lines(sqrt(VARB),rclust,lwd=3,col="black")
   lines(sqrt(VARB),rind,lwd=3, col="gray")
   lines(c(.001,10^2),c(1,1),lty=2,lwd=2,col="gray")
   
   if (!is.null(legloc)){
      legend(legloc,lwd=3,col=c("black","gray"),legend=c("by cluster","by individual"))
   }
}
par(mfrow=c(1,2))
plotSimSum(simSum[,,"MSE"],"Mean Squared Error",legloc=NULL)
plotSimSum(simSum[,,"cilen"],"Confidence Interval Length",legloc="bottomleft")

#dev.print(pdf,file="../book/graph/chStrsimRandomEffectsClusterOrdinal.pdf")

proc.time() - t0
