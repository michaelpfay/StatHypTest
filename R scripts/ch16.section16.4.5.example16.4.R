# Section 16.4.5, example 16.4
library(survival)
library(bpcp)
#library(ggplot2)
data(leuk2)
levels(leuk2$treatment)
leuk2$treatment<-factor(leuk2$treatment,levels=c("placebo","6-MP"))

cph<-coxph(Surv(time,status)~treatment, data=leuk2)


summary(cph)


leuk2$id<- factor(1:nrow(leuk2))

cph.robust<-coxph(Surv(time,status)~treatment+cluster(id), data=leuk2)

summary(cph.robust)



