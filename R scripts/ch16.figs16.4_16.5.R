distributionForms<-function(t,F, excludeTime0=TRUE){
  n<-length(t)
  if (length(F)!=n) stop("length of F must equal length of t")
  if (t[1]!=0) stop("t[1] must be zero")
  delta<- diff(t)
  f<- c(0,diff(F)/delta)
  S<- 1-F
  Sminus<- c(1,S[-n])
  h<- f/Sminus
  H<- cumsum(h*c(0,delta))
  o<- F/S
  #o<-S/F
  D<-data.frame(t=t,F=F,f=f,S=S,Sm=Sminus,h=h,H=H,o=o)
  if (excludeTime0){
    D <- D[-1,]
  }
  D
}
