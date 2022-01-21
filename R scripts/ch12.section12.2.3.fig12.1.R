# Section 12.2.3 and Fig 12.1


## the idea is to compare the difference in means (ignoring cluster) 
## to difference in means of cluster means
## with random effects difference in means
##
## I want to show that the random effects picks out the right one,
## when the varw>>varb then ignore cluster, 
## and when varw<<varb use difference in mean of cluster means
## in between cases uses the random effects model. 

t0<-proc.time()
library(lme4)
library(lmerTest)
nsim<-10^4
set.seed(104301)
NN1<-NN2<- 10
#NN1<-NN2<- 5


simCluster<-function(N1=NN1,N2=NN2,mu1=0,mu2=1, mClustSet=c(10,20,40), varw=1, varb=1){
   N<-N1+N2
   m<- sample(mClustSet,N, replace=TRUE)
   b<- rnorm(N,sd=sqrt(varb))
   B<- rep(b,m)
   CID<- factor(rep(1:N,m))
   MU<- rep(c(mu1,mu2),c(sum(m[1:N1]),sum(m[(N1+1):N])))
   G<- factor(rep(c("g1","g2"),c(sum(m[1:N1]),sum(m[(N1+1):N]))))
   Y<- MU+B+ rnorm(sum(m),sd=sqrt(varw))
   data.frame(MU,B,Y,G,CID)
}



d<-simCluster()

analyze<-function(d,trueDelta=0){
   # one obs one vote
   tout<-t.test(d$Y~d$G)
   t1<-c(estimate=unname(tout$estimate[2]-tout$estimate[1]),
         lo=-tout$conf.int[2], hi=-tout$conf.int[1],
         p.value=tout$p.value)

   # one cluster one vote
   
   #theta1Obs
   uid<-unique(d$CID)
   mu<-g<- rep(NA,length(uid))
   for (i in 1:length(uid)){
       mu[i]<- mean(d$Y[d$CID==uid[i]])
       g[i]<- as.character(d$G[d$CID==uid[i]])[1]
   }
   tout2<-t.test(mu~g)
   t2<-c(estimate=unname(tout2$estimate[2]-tout2$estimate[1]),
         lo=-tout2$conf.int[2], hi=-tout2$conf.int[1],
         p.value=tout2$p.value)
   lout<-lmer(Y~ G + (1 | CID), data= d)
   s<-summary(lout)$coefficients
   Delta<- unname(s["Gg2","Estimate"])
   seDelta<- s["Gg2","Std. Error"]
   dfDelta<- s["Gg2","df"]
   pDelta<- 2*( 1- pt(abs(Delta/seDelta),df=dfDelta) )
   loDelta<- Delta - qt(0.975,df=dfDelta)*seDelta
   hiDelta<-  Delta + qt(0.975,df=dfDelta)*seDelta
   t3<-c(estimate=Delta,
         lo=loDelta, hi=hiDelta,
         p.value=pDelta)
   #if (t1["p.value"]<=0.05) browser()
   extraStats<- function(s,trueDelta=1){
      c(s, cilen=unname(s["hi"]-s["lo"]),
        cover=unname(ifelse(s["lo"]<=trueDelta & s["hi"]>=trueDelta,1,0)),
        reject=unname(ifelse(s["p.value"]<=0.05,1,0)))
   }
   t1<-extraStats(t1,trueDelta)
   t2<-extraStats(t2,trueDelta)
   t3<-extraStats(t3,trueDelta)
   
   list(t1=t1,t2=t2,t3=t3)
}

temp<-analyze(d)
VARB<- c(.01,.1,1,10)^2
out<-array(NA,c(3,7,nsim,length(VARB)),dimnames=list(c("t.test ind","t.test clust","lmer"),
                                                       names(temp$t1),
                                                       1:nsim,
                                                     paste0("varb=",VARB)
                                                     ))
for (j in 1:length(VARB)){
   for (i in 1:nsim){
      d<-simCluster(mu2=0,varb=VARB[j])
      # trueDelta= mu2-mu1
      temp<- analyze(d,trueDelta = 0)
      out[1, ,i,j]<-temp$t1
      out[2, ,i,j]<-temp$t2
      out[3, ,i,j]<-temp$t3
   }
}
   

summarizeSim<-function(out,Delta=1){

   getmean<-function(x){ mean(x["estimate",]) }
   getMSE<- function(x){ mean((x["estimate",]-Delta)^2) }
   getCover<- function(x){ mean(x["cover",])}
   getPower<- function(x){ mean(x["reject",])}
   getCIlen<- function(x){ mean(x["cilen",])}
   sout<- array(NA,c(3,length(VARB),5),
                dimnames=list(c("t.test ind","t.test clust","lmer"),
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



plotSimSum<-function(ss,TITLE="",legloc=NULL){
   rclust<- ss[2,]/ss[3,]
   rind<- ss[1,]/ss[3,]
   plot(range(sqrt(VARB)),range(c(rclust,rind)),type="n",log="xy",xlab=expression(sigma[b]),ylab="ratio over lme",main=TITLE)
    lines(sqrt(VARB),rclust,lwd=3,col="black")
   lines(sqrt(VARB),rind,lwd=3, col="gray")
   lines(c(.001,10^2),c(1,1),lty=2,lwd=2,col="gray")
   
   if (!is.null(legloc)){
      legend(legloc,lwd=3,col=c("black","gray"),legend=c("by cluster","ignores cluster"))
   }
}


#save(simSum, file="chStrsimRandomEffectsCluster.RData")
#load("chStrsimRandomEffectsCluster.RData")
#save(simSum, file="chStrsimRandomEffectsCluster5each.RData")
#load("chStrsimRandomEffectsCluster5each.RData")
VARB<- c(.01,.1,1,10)^2

par(mfrow=c(1,2))
plotSimSum(simSum[,,"MSE"],"Mean Squared Error",legloc=NULL)
plotSimSum(simSum[,,"cilen"],"Confidence Interval Length",legloc="bottomleft")

#dev.print(pdf,file="../book/graph/chStrsimRandomEffectsCluster.pdf")
#dev.print(pdf,file="../book/graph/chStrsimRandomEffectsCluster5each.pdf")

proc.time() - t0
