## example 14.13
d<-data.frame(weightLT180=c(1,0,1,0,1,0,1,0),
          heightLT6=c(1,1,1,1,0,0,0,0),
          diabetes=c(1,1,0,0,1,1,0,0),
          y=c(34,51,536,229,1,14,39,96))

d<-data.frame(weightLT180=c(1,0,1,0),
              heightLT6=c(1,1,0,0),
              y=c(34,51,1,14),
              n=c(570,280,40,110))

d<-data.frame(weightLT180=c(1,0,1,0),
              heightLT6=c(1,1,0,0),
              y=2*c(34,51,1,14),
              n=2*c(570,280,40,110))

glmwald<-function(g){
  s<-summary(g)
  beta<-s$coef[,1]
  sebeta<-s$coef[,2]
  z<- beta/sebeta
  p.value<- 2*(1-pnorm(abs(z)))
  lower<- beta - sebeta*qnorm(.975)
  upper<- beta + sebeta*qnorm(.975)
  data.frame(OR=exp(beta),lower=exp(lower),upper=exp(upper),p.value=p.value)
}


g0<-glm(cbind(y,n-y)~heightLT6,data=d,family=binomial())
summary(g0)

glmwald(g0)

g1<-glm(cbind(y,n-y)~heightLT6+weightLT180,data=d,family=binomial())
summary(g1)

glmwald(g1)
