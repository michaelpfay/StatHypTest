# Example 10.4

compareCIs<-function(x,n){
  ciexact<-binom.test(x,n)$conf.int
  t<- x/n
  ciwald<- t + c(-1,1)*qnorm(0.975)*sqrt(t*(1-t)/n)
  
  citrans<- plogis( 
    qlogis(t) + c(-1,1)*qnorm(0.975)/sqrt(t*(1-t)*n)
    )
  ci<- matrix(c(ciexact,ciwald,citrans),nrow=3,byrow=TRUE,dimnames=list(c("exact","wald","trans"),c("lo","hi")))
  ci
}


round(compareCIs(10,100),3)
