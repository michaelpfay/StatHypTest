tscPoisTest<-function(xv, xc, nv=1,nc=1, parmtype=c("VE","ratio"), alternative=c("two.sided","less","greater"), nullValue=NULL, conf.level=0.95, force.ci2sided=TRUE){
  
  parmtype<- match.arg(parmtype)
  alternative<- match.arg(alternative)
  
  if (is.null(nullValue)){
    if (parmtype=="VE"){
      nullValue<- 0
    } else if (parmtype=="ratio"){
      nullValue<- 1
    }
  }
  

  if (parmtype=="VE"){
    VE2p<-function(v,r=nc/nv){ (1-v)/(r+1-v) }
    p0<- VE2p(nullValue)
  } else if (parmtype=="ratio"){
    ratio2p<-function(beta,r=nc/nv){ beta/(r+beta) }
    p0<- ratio2p(nullValue)
  } else stop("parmtype should be either 'VE' or 'ratio'")
  
  
  m<- xv+xc
  
  ## for VE  switch one-sided alternatives so it works for binomial parameter
  if (parmtype=="VE"){
    if (alternative=="less"){
      alternative<-"greater"
    } else if (alternative=="greater"){
      alternative<-"less"
    }
  }
  bout<- binom.test(xv,m,p=p0,conf.level=conf.level, alternative = alternative)
  ## for VE  switch alternative back for description of the output
  if (parmtype=="VE"){
    if (alternative=="less"){
      alternative<-"greater"
    } else if (alternative=="greater"){
      alternative<-"less"
    }
  }

  
  # often we want a one-side p-value, but we want the CI to be two-sided
  if (force.ci2sided){
    ci<- binom.test(xv,m,p=p0,conf.level=conf.level, alternative = "two.sided")$conf.int
  } else {
    ci<- bout$conf.int    
  }

  Method<-c("Two-Sample Conditional Poisson Test")
  if (force.ci2sided & alternative!="two.sided"){
    Note<- " (CI two-sided, but p-value one-sided)"
    Method<- paste0(Method, Note)
  } 
  
  if (parmtype=="VE"){
    p2VE<-function(p,r=nc/nv){ 1 - ( p*r/(1-p) ) }
    cint<- p2VE(c(ci[2],ci[1]))
    estimate<- p2VE(xv/m)
    names(estimate)<-"VE"
    names(nullValue)<-"VE null"
  } else if (parmtype=="ratio"){
    p2ratio<-function(p,r=nc/nv){ p*r/(1-p)  }
    cint<- p2ratio(c(ci[1],ci[2]))
    estimate<- p2ratio(xv/m)
    names(estimate)<-"ratio"
    names(nullValue)<-"ratio null"
  } else stop("parmtype should be either 'VE' or 'ratio'")

  attr(cint,"conf.level")<- conf.level

  r<- nc/nv
  names(r)<-"nc/nv"
    dname <- paste("xv=",deparse1(substitute(xv)), " and xc=", 
                   deparse1(substitute(xc)))
    
  out <- list(statistic = c(xv=xv,xc=xc), parameter = r, p.value = bout$p.value, 
               conf.int = cint, estimate = estimate, null.value = nullValue, 
               alternative = alternative, method = Method, 
               data.name = dname)
  class(out) <- "htest"
  out
}

# See Sadoff, et al, 2021,  Table~2 for >=14 days after administration 
tscPoisTest(116,348,nv=3116.6,nc=3096.1,nullValue = 0.30, alternative="greater")
# compare confidence interval to poisson.exact version
library(exactci)
1- poisson.exact(c(116,348),c(3116.6,3096.1))$conf.int[2:1]