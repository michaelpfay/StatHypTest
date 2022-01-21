x1<- 10
n1<-63
x2<-67
n2<-69

x1/n1
x2/n2
x2/n2 - x1/n1


(x2/n2)/(x1/n1)

((n1-x1)/n1)/((n2-x2)/n2)


dbinom(0,10,.01)*dbinom(0,10,.0001)

library(exact2x2)

########################
# FAKE DATA
###########################

# Example 7.3
x1<-2
n1<-11
x2<-6
n2<-8
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="central")
p1<-x1/n1
p2<-x2/n2
p2*(1-p1)/(p1*(1-p2))

fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="central")
minlike<-fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="minlike")
minlike
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="minlike", conf.level=1-minlike$p.value)

binomMeld.test(x1,n1,x2,n2, parmtype="d")
binomMeld.test(x1,n1,x2,n2, parmtype="r")
# typo in the book, it should have an upper limit of 382.95
binomMeld.test(x1,n1,x2,n2, parmtype="o")

plotHypergeometric<-function(x1,n1,x2,n2){

xmin<- max(0,(x1+x2)-n1)
xmax<- min(x1+x2,n2)
X2<-xmin:xmax
f<-dhyper(X2,n2,n1,x1+x2)
sum(f)
## check p-values
## Fisher-Irwin
sum(f[f<=f[X2==x2]])
# Central FET
2*sum(f[X2>=x2])
plot(c(xmin-1,xmax+1),range(f),type="n",xlab=expression(X[2]),ylab="probability",axes=FALSE)
axis(1,X2)
axis(2,las=2)
box()
drawfi<-function(xi,fi){
  if (xi==x2){
    COL<-"gray"
  } else if (fi<=f[X2==x2]){
    COL<-"black"
  } else COL<-"white"
  polygon(c(xi-.5,xi-.5,xi+.5,xi+.5,xi-.5),
          c(0,fi,fi,0,0),col=COL,border="black")
}
for (i in 1:length(X2)){
  drawfi(X2[i],f[i])
}
}
# Fig 7.1
plotHypergeometric(x1,n1,x2,n2)

#dev.print(pdf,file="../graph/ch2binQuickCals.pdf")



###################################
# Central FET p < Fisher-Irwin
###################################
x1<-0
n1<-10
x2<- 22
n2<- 100

plotHypergeometric(0,10,22,100)


fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="central")$p.value
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),tsmethod="minlike")$p.value


#######################################
# Non-unified example
#######################################
# Example 2.5
x1<-7
n1<-262
x2<-30
n2<-494
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1,tsmethod="central")$p.value
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1-.01,tsmethod="central")$p.value
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1.01,tsmethod="central")$p.value

fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1,tsmethod="minlike")$p.value
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1.01,tsmethod="minlike")$p.value
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1-.01,tsmethod="minlike")$p.value


######################################
# PREVAIL II
#####################################
# Example 7.4
x1<-0
n1<-24
x2<-4
n2<-20
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1,tsmethod="central")
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1,tsmethod="minlike")
fisher.exact(matrix(c(x2,n2-x2,x1,n1-x1),2,2),or=1,alternative = "greater")


