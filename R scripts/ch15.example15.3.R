## Fake  Measles Vaccine  Example 15.3


getTables<-function(VE=.80,prevlo=10^(-5), prevhi=10^(-3), propHi=.2,vlo=.70, vhi=.95){
  p<-array(1:8,dim=c(2,2,2),dimnames=list(
    c("z=0","z=1"),c("Y=0","Y=1"),c("x=0","x=1")
  ))
  
  

  
  
  ## x=0 low
  p0.0<- (1-vlo)*(1-propHi)
  p1.0<-  vlo*(1-propHi)
  p010<- prevlo*p0.0
  p110<- p1.0*(1-VE)*prevlo

  #p(0,1,0)<- p010
  #p(0,0,0)<- p0.0 - p010
  #p(1,1,0)<-p110
  #p(1,0,0)<-p1.0 - p110
  
  p[1,2,1]<- p010
  p[1,1,1]<- p0.0 - p010
  p[2,2,1]<-p110
  p[2,1,1]<-p1.0 - p110
  
  
  ## x=1 hi
  p0.1<- (1-vhi)*(propHi)
  p1.1<-  vhi*(propHi)
  p011<- prevhi*p0.1
  p111<- p1.1*(1-VE)*prevhi

  #p(0,1,1)<- p011
  #p(0,0,1)<- p0.1 - p011
  #p(1,1,1)<-p111
  #p(1,0,1)<-p1.1 - p111
  
  p[1,2,2]<- p011
  p[1,1,2]<- p0.1 - p011
  p[2,2,2]<-p111
  p[2,1,2]<-p1.1 - p111
  
  p
}

p<-getTables()


(p[2,2,1]/(p[2,1,1]+p[2,2,1]))/ 
(p[1,2,1] /(p[1,1,1]+p[1,2,1]) )


(p[2,2,2]/(p[2,1,2]+p[2,2,2]))/ 
  (p[1,2,2] /(p[1,1,2]+p[1,2,2]) )

pp<-apply(p,c(1,2),sum)

(pp[2,2]/(pp[2,1]+pp[2,2]))/ 
  (pp[1,2] /(pp[1,1]+pp[1,2]) )

n<- 10^8 * p
n


tot<-apply(n,c(1,3),sum)
tot


