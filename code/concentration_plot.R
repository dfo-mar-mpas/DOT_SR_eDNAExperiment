#plot up eDNA concentration vs duration in RNAlater in the DOT samplers
library(dplyr)
library(ggplot2)

#Read in the csv 
conc <- read.csv("data/eDNA_Concentrations.csv", header = T) %>% group_by(Method)

ggplot(conc, aes(x=DurationPreserved, y=as.numeric(Concentration), color=Method))+
  geom_point(size=3)+
  geom_smooth(se = F)+
  ylab(label = "DNA concentration ng/ul")+
  xlab(label = "Days preserved")+
  expand_limits(x=0,y=0)+
  #scale_x_continuous(expand=c(0,0), limits=c(0,NA))+
  scale_y_continuous(expand=c(0,0), limits=c(0,NA))+
  theme_bw()+
  theme(text=element_text(size=18))


ggsave(filename = "DOTvsSR_preservation.png",plot = last_plot(),device = "png",
       path = "figures/",width = 12, height = 8, units = "in", dpi = 400, bg = "white")
