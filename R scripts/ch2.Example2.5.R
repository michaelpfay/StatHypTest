library(exact2x2)

x<- matrix(c(7,262-7,30,494-30),2,2)
fisher.exact(x,or=1,tsmethod="central")
fisher.exact(x,or=1,tsmethod="minlike")




fisher.test(x,or=.99)

fisher.test(x,or=1)


fisher.test(x,or=1.01)



exact2x2Plot(matrix(c(7,262-7,30,494-30),2,2),orRange=c(0.175,1.02),tsmethod="minlike")
exact2x2Plot(matrix(c(7,262-7,30,494-30),2,2),orRange=c(0.17,0.190),tsmethod="minlike")
lines(rep(0.1773,2),c(0,.05),col="red")

exact2x2Plot(matrix(c(7,262-7,30,494-30),2,2),orRange=c(0.9,1.2),tsmethod="minlike")

exact2x2Plot(matrix(c(7,262-7,30,494-30),2,2),orRange=c(0.95,1.05),tsmethod="minlike")
exact2x2Plot(matrix(c(7,262-7,30,494-30),2,2),orRange=c(0.95,1.05),tsmethod="minlike",ylim=c(0.0495,0.0505))

lines(rep(0.993,2),c(0,.05),col="blue")
lines(rep(1.006,2),c(0,.05),col="brown")
lines(rep(1.014,2),c(0,.05),col="green")