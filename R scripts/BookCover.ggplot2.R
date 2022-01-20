library(bpcp)
data(leuk2)
leuk2fit <- bpcpfit(Surv(time, status)~treatment, data=leuk2)
leuk2tidy <- tidykmciLR(leuk2fit)

fill.colors<- c(gray(.1),gray(.5))
line.colors<- c(gray(.1),gray(.5))

ggplot(leuk2tidy, 
  aes(x = time, y = surv, ymin = lower, ymax = upper, col = group)) + 
  geom_line(show.legend=FALSE,lwd=1,aes(linetype=group)) + 
  geom_ribbon(alpha = .2, aes(fill=group,linetype=group)) + 
  xlab(" ") + ylab(" ") + 
  ggtitle("") +   theme_void() +
  theme(legend.position="None") +
  scale_fill_manual(values=fill.colors)+
  scale_color_manual(values=line.colors)+
  scale_linetype_manual(values=c("solid","dashed"))

#ggsave(filename="cover900.tiff",dpi=900)
#ggsave(filename="cover300.tiff",dpi=300)
