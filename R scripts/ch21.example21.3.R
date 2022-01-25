# example 21.3
## Bayes Factor example with Large Binomial data

x<-52263471
n<-104490000

## Pr[X=x | theta=.5]
dbinom(x,n,.5)
exp(lchoose(n,x)+n*log(.5))

## Pr[ X=x | prior=Beta(a,b) ] 
## Use beta distribution to know that
##      int_0^1 t^(x-1) (1-t)^(y-1) dt = gamma(x)*gamma(y)/gamma(x+y) 
PrXgivenPrior<-function(a,b){
    exp(lchoose(n,x)+lgamma(x+a)+lgamma(n-x+b) - lgamma(n+a+b) )
}

dbinom(x,n,.5)/PrXgivenPrior(1,1)

dbinom(x,n,.5)/PrXgivenPrior(.5,.5)


