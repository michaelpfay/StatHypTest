negbinom.vs.binom<-function(x,n){
   # gives one-sided p-values
   c(p.binom=pbinom(x,n,.5),
     p.nbinom=1-pnbinom(n-x-1,x,.5))
}

#negbinom.vs.binom(3,12)
negbinom.vs.binom(6,20)