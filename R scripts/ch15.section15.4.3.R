# Section 15.4.3, Table 15.1


n<-function(a,b){
  if (a==0 & b==0) out<- 2385
  if (a==0 & b==1) out<- 9663
  if (a==1 & b==0) out<- 34
  if (a==1 & b==1) out<- 12
  out
}

m0. <- 11514
m1. <- 74
M<- m0.  + m1.
N<- n(0,0)+n(0,1)+n(1,0)+n(1,1)

m01<- m0. - (M/N)*n(0,0)
m11<- m1. -(M/N)*n(1,0)


Rtilde1<-  (n(1,1)/(n(0,1)+n(1,1)) )/(m11/(m01+m11))
Rtilde1
## eqn 3: there is a typo in Sommer and Zeger
##
Rtilde3<- n(1,1)*(m1.*(M/N)-n(1,0))
Rtilde3
## I think it should be...
Rtilde4<- n(1,1)/ ((N/M)*m1. - n(1,0))
Rtilde4


### 
pi<- 9675/12094
theta00 <- 74/11588
theta10.T1eq0<-  34/2419
theta00.T1eq1 <-  (theta00 - (1-pi)* theta10.T1eq0 )/(pi)
theta00.T1eq1
theta11.T1eq1<- 12/9675
theta11.T1eq1



theta11.T1eq1/theta00.T1eq1

0.00124/0.00447

theta11.T1eq1-theta00.T1eq1