# Section 18.1
# Sequential t-test Example

t.test.fast<-function(x){
  n<-length(x)
  stderr<- sqrt(var(x)/n)
  2*(1-pt(abs(mean(x)/stderr),df=n-1))
}


simulateRepeat.ttest <- function(Nmax,
                              nsim = 1e3,
                              alpha.two.sided = 0.05) {
  N <- 1:Nmax
  reject <- rejectBonf<-rep(NA, nsim)
  for (i in 1:nsim) {
    Y <- rnorm(Nmax)
    p<-rep(1,Nmax)
    for (j in 2:Nmax){
      p[j]<- t.test.fast(Y[1:j])
    }
    reject[i] <- any(p <= alpha.two.sided)
    rejectBonf[i] <- any(p <= (alpha.two.sided/(Nmax-1)) )
  }
  list(power=sum(reject) / nsim, powerBonf=sum(rejectBonf)/nsim, minReject=2) 
}

set.seed(33394928)
simulateRepeat.ttest(101,nsim=1e4,alpha.two.sided=0.05)

set.seed(9694)
simulateRepeat.ttest(1001,nsim=1e4,alpha.two.sided=0.05)

