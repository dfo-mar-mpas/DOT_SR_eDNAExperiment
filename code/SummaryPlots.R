## Plotting DOT comparison data -----

#load libraries ---
library(tidyverse)

#functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#load data --- 
dot_data_raw <- read.csv("data/12S_Taxonomy_Table_FilteredFinal_grouped.csv")

#format to long 
dot_df <- dot_data_raw%>%
          gather(key = "Sampler",val="reads",3:48)%>%
          mutate(type=ifelse(grepl("DOT",Sampler),"Autonomous","FAS"),
                 week=substrRight(Sampler,1),
                 rep=case_when(grepl("1003",Sampler)~1,
                               grepl("1004",Sampler)~2,
                               grepl("1005",Sampler)~3,
                               grepl("SR1",Sampler)~1,
                               grepl("SR2",Sampler)~2,
                               grepl("SR3",Sampler)~3,
                               TRUE ~ NA)) #only a catch if it didn't work: sum(is.na(dot_df$rep)) == 0
                             

dot_df_week_rich <- dot_df%>%
                    filter(reads>0)%>%
                    group_by(type,week,rep)%>%
                    summarize(nspecies=length(unique(Sp)))%>%
                    ungroup()%>%
                    group_by(type,week)%>%
                    summarize(mn=mean(nspecies),
                              sd=sd(nspecies))%>%
                    ungroup()%>%
                    mutate(se=sd/sqrt(3),
                           se=ifelse(week>3 & type == "SR",sd,se))%>%
                    data.frame()
  

p1 <- ggplot(data=dot_df_week_rich,aes(x=week,y=mn))+
  geom_errorbar(aes(ymin=mn-se,ymax=mn+se,group=type,col=type),
                linetype = "dotted",
                position=position_dodge(width=0.75),
                width=0.2,show.legend = FALSE)+
  geom_point(position=position_dodge(width=0.75),size=3,aes(group=type,fill=type),pch=21,col="black")+
  geom_text(data=data.frame(week = c(1,9),
                            mn = rep(range(dot_df_week_rich$mn)[1]-2,2), 
                            label = c("Aug 11 2023","Oct 16 2023")),aes(label=label))+
  theme_bw()+
  labs(x="Sample week",y="Species detected",fill="")+
  theme(legend.position = c(0, 1),  # Position legend in the top left
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.title = element_blank())

ggsave("figures/AverageSpeciesRichness.png",p1,width=7.5,height=4.5,units="in",dpi=300)


dot_df_week_rich_total <- dot_df%>%
  filter(reads>0)%>%
  group_by(type,week)%>%
  summarize(nspecies=length(unique(Sp)))%>%
  ungroup()%>%
  data.frame()


p2 <- ggplot(data=dot_df_week_rich_total,aes(x=week,y=nspecies))+
  geom_path(lwd=0.25,aes(group=week),lty=2,show.legend = FALSE)+
  geom_point(size=2,aes(group=type,fill=type),pch=21,col="black")+
  geom_text(data=data.frame(week = c(1,9),
                            nspecies = rep(range(dot_df_week_rich_total$nspecies)[1]-2,2), 
                            label = c("Aug 11 2023","Oct 16 2023")),aes(label=label))+
  theme_bw()+
  labs(x="Sample week",y="Species detected",fill="")+
  theme(legend.position = c(0, 1),  # Position legend in the top left
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.title = element_blank())

ggsave("figures/TotalRichness_Comparison.png",p2,width=7.5,height=4.5,units="in",dpi=300)
