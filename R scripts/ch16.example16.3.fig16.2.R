# example 16.3, fig 16.2
library(bpcp)
library(ggplot2)
data(leuk2)
timePlacebo<- leuk2[leuk2$treatment=="placebo","time"]
statusPlacebo<- leuk2[leuk2$treatment=="placebo","status"]

time6mp<- leuk2[leuk2$treatment=="6-MP","time"]
status6mp<- leuk2[leuk2$treatment=="6-MP","status"]


library(survival)

sfit<-survfit(Surv(timePlacebo,statusPlacebo)~1)
str(sfit)


fitp<-bpcp(timePlacebo,statusPlacebo)
fit6mp<-bpcp(time6mp,status6mp)

plot(fit6mp)
lines(fitp, lwd=3)
legend("topright",legend=c("6-MP","placebo"),lty=c(1,1),lwd=c(1,3))

#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenDiffSt.pdf")



bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,20,parmtype="ratio")

bpcp2samp(leuk2$time,leuk2$status,leuk2$treatment,20,parmtype="difference")



### Prettier!!

library(bpcp)
library(ggplot2)
library(tidyr)

data(leuk2)
timePlacebo<- leuk2[leuk2$treatment=="placebo","time"]
statusPlacebo<- leuk2[leuk2$treatment=="placebo","status"]

time6mp<- leuk2[leuk2$treatment=="6-MP","time"]
status6mp<- leuk2[leuk2$treatment=="6-MP","status"]

fitp<-bpcp(timePlacebo,statusPlacebo)
fit6mp<-bpcp(time6mp,status6mp)

#Function to "tidy" output from bpcp function (kmciLR object)
#Input:
#kmciLR.data: Vector of kmciLR objects - output from previously fit bpcp functions
#labels: group (treatment) variables for the different kmciLR objects. Default is NA. Order should be the same as
#the order of the kmciLR objects. Vector of strings preferred
#output is a dataframe with time, survival, upper and lower CI limits, and group 
tidy.kmciLR <- function(kmciLR.data, labels = NA) {
  
  #Create empty dataframe to store output
  tidyout <- data.frame(NULL)
  #Determine number of groups
  num <- (length(kmciLR.data)/12)
  
  #If no labels are provided and there is more than one kmciLR object, create numeric labels
  if (num > 1 && is.na(labels)) {
    labels <- c(1:num)
  }
  
  #Error handling: Check that number of labels provided equals number of groups
  else if (num > 1 && (length(labels) != num)) {
    stop("Number of labels must equal number of kmciLR objects")
  }
  
  #For each group, "tidy" the output into a new dataframe
  for (i in 1:num) {
    new <- with(kmciLR.data[((i*12)-11):(i*12)], data.frame(time = sort(c(L, R)), 
                                                            surv = rep(surv, each = 2), 
                                                            lower = rep(lower, each = 2), 
                                                            upper = rep(upper, each = 2)))
    #Add group variable 
    if (!all(is.na(labels))) {
      new$Group <- labels[i]
    }
    #Add to tidy dataframe
    tidyout <- rbind(tidyout, new)
  }
  #Factors the group variable to account for numbers being entered instead of strings
  if (!all(is.na(labels))) {
    tidyout$Group <- as.factor(tidyout$Group)
  }
  return(tidyout)
  
}

#Example 1: leuk2 data, 2 groups, labels provided
fit_data <- c(fitp, fit6mp)
tidyout <- tidy.kmciLR(fit_data, labels = c("Placebo", "6-MP"))
ggplot(tidyout, aes(x = time, y = surv, ymin = lower, ymax = upper, linetype = Group, color=Group, fill=Group)) + 
  geom_line() + geom_ribbon(alpha = .2, colour=gray(.5))

ggplot(tidyout, aes(x = time, y = surv, ymin = lower, ymax = upper, linetype = Group)) + 
  geom_line() + geom_ribbon(alpha = .2, colour=gray(.5))

#dev.print(pdf,file="S:/BRB/Staff Folders/Mike Fay/book/ASHT/book/graph/chCenDiffSt.pdf")


#Example 2: one group, with and without labels
#Shows modifications to plot
ggplot(tidy.kmciLR(fitp), aes(x = time, y = surv, ymin = lower, ymax = upper)) + 
  geom_line() + geom_ribbon(alpha = .2) + ggtitle ("KM-plot for Placebo") +xlab("Time (Weeks)") + ylab("Survival")

ggplot(tidy.kmciLR(fitp, "Placebo"), aes(x = time, y = surv, ymin = lower, ymax = upper, linetype = Group)) + 
  geom_line() + geom_ribbon(alpha = .2) + ggtitle ("KM-plot for Placebo") +xlab("Time (Weeks)") + ylab("Survival")

