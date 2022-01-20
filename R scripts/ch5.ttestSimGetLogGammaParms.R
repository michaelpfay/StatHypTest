library(VGAM)
SHAPE<-.2
rootfunc<-function(s){
    # we want scale parameter so that Pr[ -1 < X < 1]=.5
    loc<-1-qlgamma(.75,shape=SHAPE,scale=s)
    plgamma(1,location=loc,shape=SHAPE,scale=s) - 
       plgamma(-1,location=loc,shape=SHAPE,scale=s) - 0.5
}

rootfunc(.1)
rootfunc(200)
s<-uniroot(rootfunc,c(.1,200))$root
rootfunc(s)
loc<- 1-qlgamma(.75,shape=SHAPE,scale=s)
# checks
qlgamma(.75,location=loc,shape=SHAPE,scale=s)
qlgamma(.25,location=loc,shape=SHAPE,scale=s)


plgamma(-1,location=loc,shape=SHAPE,scale=s)
plgamma(1,location=loc,shape=SHAPE,scale=s)


### ERICA and I spoke about shifting the mean to equal 0
### but we decided against it
###    here is shifting the mean to 0....
### Now shift location so that E(Y)=0
### E(Y)  = loc + scale*digamma(shape)
#EYold<-loc+s*digamma(SHAPE)
#loc<- loc -EYold
#EY<-loc+s*digamma(SHAPE)
#################################################

logGammaDist<-function(x){ 
    dlgamma(x,location=loc,shape=SHAPE,scale=s)
}
logGammaMean<-loc+s*digamma(SHAPE)
