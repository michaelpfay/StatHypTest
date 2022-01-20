# Ch 6 example 6.4
# McNemar's test example
B<- 9
C<- 2



ci<-binom.test(B,B+C)$conf.int
binom.test(B,B+C)

tr<-function(p){ p/(1-p) }

p<- B/(B+C)

tr(p)

tr(ci)
