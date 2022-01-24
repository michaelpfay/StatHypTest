# Tables 18.1 and 18.2, example 18.1
library(gsDesign)


createTablePocock<-function(overallAlpha=0.05,K=c(1,2,3,5,10,20,30)){
  Abonf<- Apocock<-Cpocock<-rep(NA,length(K))
  ALPHA<-overallAlpha
  for (i in 1:length(K)){
    if (K[i]==1){
      Cpocock[i]<- qnorm(1-ALPHA/2)
      Apocock[i]<-  2*(1-pnorm(Cpocock[i]))
    } else if (K[i]<=30){
      g<-gsDesign(k=K[i],test.type=2, alpha=ALPHA/2,sfu="Pocock")
      Cpocock[i]<-g$upper$bound[1]
      Apocock[i]<- 2*(1-pnorm(Cpocock[i]))
    } else stop("K must be between 2 and 30")
  }
  Abonf<- ALPHA/K
  data.frame(K,Abonf,Apocock,Cpocock)
}



d<-createTablePocock()
dimnames(d)[[1]]<-NULL
library(xtable)
print(xtable(d,digits=c(0,0,4,4,2)),include.rownames=FALSE)


#g3<-gsDesign(k=3,test.type=2,sfu="Pocock")
# Calculation for example
g5<-gsDesign(k=5,test.type=2,sfu="Pocock")
g5


probStop<-g5$upper$spend + g5$lower$spend
cumAlpha<-cumsum(g5$upper$spend + g5$lower$spend)
names(probStop)<-names(cumAlpha)<- 1:5
round(probStop*100,1)
cumAlpha
# with power .90
g5$upper$prob[,2]
# check 
sum(g5$upper$prob[,2])



creatTableOF<-function(overallAlpha=0.05,K=c(2,3,4,5)){
  cOF<-rep(NA,length(K))
  ALPHA<-overallAlpha
  nomAlpha<-matrix(NA,length(K),max(K))
  for (i in 1:length(K)){
      g<-gsDesign(k=K[i],test.type=2, alpha=ALPHA/2,sfu="OF")
      cOF[i]<- g$upper$bound[K[i]]
      t<-  1:K[i]/K[i]
      Zbound<- cOF[i]/sqrt(t)
      nomAlpha[i,1:K[i]]<- 2*(1-pnorm(Zbound))
  }

  data.frame(K,cOF,nomAlpha)
}

dof<-creatTableOF()
print(xtable(dof,digits=c(0,0,2,6,3,3,3,3)),include.rownames=FALSE)

# Calculation for example
g5<-gsDesign(k=5,test.type=2,sfu="OF")
g5


probStop<-g5$upper$spend + g5$lower$spend
cumAlpha<-cumsum(g5$upper$spend + g5$lower$spend)
names(probStop)<-names(cumAlpha)<- 1:5
round(probStop*100,3)
cumAlpha
# with power .90
round(100*g5$upper$prob[,2],1)
# check 
sum(g5$upper$prob[,2])


#g30<-gsDesign(k=30,test.type=2,sfu="OF")

