# Behrens-Fisher distribution

pbf<-function(q,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8)){
    if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
        R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
    } else if (is.null(s1) & is.null(s2)){
        # do nothing R<-R
    } else stop("supply one of R or (s1 and s2)")

    sinR<-sin(R)
    cosR<-cos(R)
    
    # if sinR=0 then cosR=1 and D=T2
    # and if cosR=0 then sinR=1 and D=T1
    if (abs(sinR)<epsilon){
        out<-pt(q,n2-1)
    } else if (abs(cosR)<epsilon){
        out<-pt(q,n1-1)
    } else if (cosR>0 & cosR<1 & sinR>0 & sinR<1){
    ## usual case       
        ifunc<-function(u,Y){
            pt((Y+u*sinR)/cosR,n2-1)*dt(u,n1-1)
        }
        nq<-length(q)
        out<-rep(NA,nq)
        for (i in 1:nq){
            out[i]<-integrate(ifunc,-Inf,Inf,Y=q[i])$value
        }
     }

     out
}
#pbf(c(.9,.95,.975),12,6,R=pi/2-0.00001)
#pt(c(.9,.95,.975),11)
#z<-pbf(c(.345,.45,.6,7),12,6,s1=10,s2=1)
#z


qbf<-function(p,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8)){
    if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
        R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
    } else if (is.null(s1) & is.null(s2)){
        # do nothing R<-R
    } else stop("supply one of R or (s1 and s2)")

    ## if R=0 then cos(R)=1 and D ~ T2
    if (abs(R)<epsilon){
        out<-qt(p,n2-1)
    } else if (abs(R-pi/2)<epsilon){
        out<-qt(p,n1-1)
    } else {
        rootfunc<-function(q,P){
            P-pbf(q,n1,n2,R)
        }
        np<-length(p)
        out<-rep(NA,np)
        INTERVAL<-c(qt(epsilon,1),-qt(epsilon,1))
        for (i in 1:np){
            out[i]<-uniroot(rootfunc,interval=INTERVAL,P=p[i])$root
        }
        out
    }
    out
}

#pbf(2.91999,2,2,R=0)
#qbf(c(.95),5,8,R=pi/4)
#qbf(c(.05),4,15,R=pi/8)

