# example 13.8

## Simulation of correlated effects

## Create data 
#library(mvtnorm);set.seed(1);n<-200;k<-8;rho<-0.9
#R<-diag(k);R[R==0]<- rho;X<- rmvnorm(n,sigma=R);Y<-X[,1] + rnorm(n);d<-data.frame(X,Y)

library(multcomp)
simMVnorm <- function(n, k, rho, seed) {
  set.seed(seed)
  #n<- 200
  #k<-8
  #rho<-0.9
  R <- diag(k)
  R[R == 0] <- rho
  X <- rmvnorm(n, sigma = R)
  Y <- X[, 1] + rnorm(n)
  d <- data.frame(X, Y)
  xvars <- paste0("X", 1:k)
  fmla <- as.formula(paste("Y ~ ", paste(xvars, collapse = "+")))
  lout <- lm(fmla, data = d)
  slout <- summary(lout)
  plinear <- slout$coef[-1, "Pr(>|t|)"]
  C <- diag(k + 1)[-1, ]
  sout <- summary(glht(lout, linfct = C))
  pmaxt <- sout$test$pvalues
  
  pj <- rep(NA, k)
  for (j in 1:k) {
    fmj <- as.formula(paste("Y~", paste0("X", j)))
    pj[j] <- summary(lm(fmj, data = d))$coef[2, "Pr(>|t|)"]
  }
  pj
  pholm <- p.adjust(pj, "holm")
  
  
  ###########################
  # Repeat without x1
  
  
  
  xvars <- paste0("X", 2:k)
  fmla <- as.formula(paste("Y ~ ", paste(xvars, collapse = "+")))
  lout <- lm(fmla, data = d)
  slout <- summary(lout)
  plinearM1 <- c(NA, slout$coef[-1, "Pr(>|t|)"])
  C <- diag(k)[-1, ]
  sout <- summary(glht(lout, linfct = C))
  pmaxtM1 <- c(NA, sout$test$pvalues)
  pmaxtM1
  
  
  tab <- matrix(c(plinear, pmaxt, pj, pholm, plinearM1, pmaxtM1),
                k,
                6,
                dimnames = list(
                  paste0("X", 1:k),
                  c("lmWald", "lmMaxt", "pj", "pholm", "lmWaldM1", "lmMaxtM1")
                ))
  tab
}

simMVnorm(200,8,.9,1)
stab0<-simMVnorm(200,8,0,1)
round(stab0,4)