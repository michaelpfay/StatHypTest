library(exact2x2)

fisher.exact(matrix(c(7,262-7,30,494-30),2,2),or=1,tsmethod="central")
fisher.exact(matrix(c(7,262-7,30,494-30),2,2),or=1,tsmethod="minlike")



fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=0.177)



fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=0.178)

fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=.99)

fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=1)


fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=1.01)

fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=1.01)
fisher.test(matrix(c(7,262-7,30,494-30),2,2),or=1.013)


binom.exact(2,12,tsmethod="minlike",plot=TRUE,xlim=c(.4,.6))
exactbinomPlot(2,12,tsmethod="minlike",xlim=c(.4,.6),ylim=c(.02,.08))