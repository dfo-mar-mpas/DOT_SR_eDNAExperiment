## Plotting DOT comparison data -----

#load libraries ---
library(tidyverse)
library(patchwork)

#functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#load data --- 
dot_data_raw <- read.csv("data/12S_Revised_ASVsGrouped.csv") #fish
dot_data_raw_invert <- read.csv("data/DOTCOI_feature_table_export_filtered.csv")

#format to long 
dot_df <- dot_data_raw%>%
          gather(key = "Sampler",val="reads",3:49)%>%
          mutate(type=ifelse(grepl("X10",Sampler),"Autonomous","FAS"),
                 type=ifelse(grepl("Control",Sampler),"FAS",type),
                 week=ifelse(grepl("X10",Sampler),substrRight(Sampler,1),substr(Sampler,3,3)),
                 rep=case_when(grepl("1003",Sampler)~1,
                               grepl("1004",Sampler)~2,
                               grepl("1005",Sampler)~3,
                               grepl("Control",Sampler)~1,
                               TRUE ~ NA),
                 marker="12S")%>% #only a catch if it didn't work: sum(is.na(dot_df$rep)) == 0
          rename(confed = PercentMatch,
                 Sp = BLAST)

dot_df_invert <- dot_data_raw_invert%>%
                 gather(key = "Sampler",val="reads",7:53)%>%
                 mutate(type=ifelse(grepl("DOT",Sampler),"Autonomous","FAS"),
                         week=substrRight(Sampler,1),
                         rep=case_when(grepl("1003",Sampler)~1,
                                       grepl("1004",Sampler)~2,
                                       grepl("1005",Sampler)~3,
                                       grepl("SR1",Sampler)~1,
                                       grepl("SR3",Sampler)~2,
                                       grepl("SR4",Sampler)~3,
                                       grepl("SR5",Sampler)~4,
                                       TRUE ~ NA),
                        marker="CO1")%>%
                  rename(Sp=Species,confed = 6)

#calculate weekly richness

## Data prep --------------

rich_df <- rbind(dot_df,dot_df_invert%>%select(names(dot_df)))%>% #number of speci
            filter(reads>0)

dot_df_rich_rep <- rich_df%>% #per replicate richness
                  group_by(marker,type,week,rep)%>%
                  summarize(nspecies=length(unique(Sp)))%>%
                  ungroup()%>%
                  data.frame()

dot_df_week_rich <- dot_df_rich_rep%>% #average replicate richness
                    group_by(marker,type,week)%>%
                    summarize(mn=mean(nspecies),
                              sd=sd(nspecies))%>%
                    ungroup()%>%
                    mutate(se=sd/sqrt(3),
                           se=ifelse(week>3 & type == "SR",sd,se))%>%
                    data.frame()%>%
                    mutate(rep=NA)

dot_df_week_rich_overlap <- rich_df%>% #overlap amongst species detected per week. 
                      group_by(marker,type,week)%>%
                      distinct(Sp,.keep_all=TRUE)%>%
                      ungroup()%>%
                      group_by(marker,week,Sp)%>%
                      summarize(count=n())%>%
                      ungroup()%>%
                      group_by(marker,week)%>%
                      summarize(shared = sum(count==2),
                                count=n(),
                                overlap = shared/count)%>%
                      ungroup()%>%
                      data.frame()%>%
                      mutate(rep=NA,se=NA,sd=NA)

#Assemble plot data
plotvals <- c("marker","type","week","rep","xval","sd","se")

plotdata <- rbind(dot_df_rich_rep%>%
                    rename(xval=nspecies)%>%
                    mutate(sd=NA,se=NA,type=paste(type,"replicate",sep="_"))%>%
                    select(plotvals),
                  
                  dot_df_week_rich%>%
                    rename(xval=mn)%>%
                    dplyr::select(plotvals),
                  
                  dot_df_week_rich_overlap%>%
                    rename(Shared = shared,
                           Total=count)%>%
                    gather(key = "type", value = "xval", Shared, Total)%>%
                    mutate(sd=NA,se=NA)%>%
                    select(plotvals)
                  
)%>%mutate(type2=gsub("_replicate","",type))


#12S marker plot --------------
p1_a <- ggplot()+
  #per replicate points
  geom_point(data=plotdata%>%
               filter(marker=="12S",type%in%c("Autonomous_replicate","FAS_replicate")),
             aes(x=week,y=xval,col=type2),size=0.5, position=position_dodge(width=0.75),show.legend = FALSE)+
  
  #mean among replicates
  geom_errorbar(data=plotdata%>%
                  filter(marker=="12S",type %in% c("Autonomous","FAS")),
                aes(x=week,y=xval,ymin=xval-se,ymax=xval+se,group=type2,col=type2),
                #linetype = "dotted",
                position=position_dodge(width=0.75),
                width=0.2,show.legend = FALSE)+
  geom_point(data=plotdata%>%filter(marker=="12S",type %in% c("Autonomous","FAS")),
             position=position_dodge(width=0.75),size=3,
             aes(x=week,y=xval,group=type2,fill=type2),pch=21,col="black")+
  
  #shared species points and line
  geom_line(data=plotdata%>%
              filter(marker=="12S",type=="Shared")%>%
              mutate(marker="12S"),aes(x=week,y=xval,group=1),lty=2,col="grey60",lwd=0.6)+
  geom_point(data=plotdata%>%filter(marker=="12S",type=="Shared"),
             aes(x=week,y=xval),shape=21,fill="grey60",size=1.5)+
  
  #total species 
  geom_line(data=plotdata%>%
              filter(marker=="12S",type=="Total")%>%
              mutate(marker="12S"),aes(x=week,y=xval,group=1),lty=2,col="grey30",lwd=0.6)+
  geom_point(data=plotdata%>%filter(marker=="12S",type=="Total"),
             aes(x=week,y=xval),shape=21,fill="grey80",size=1.5)+
  
  facet_grid(marker~.)+
  theme_bw()+
  labs(x="Sample week",y="Species detected",fill="")+
  theme(legend.position = c(0, 1),  # Position legend in the top left
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,26))

#CO1 marker plot --------------
p1_b <- ggplot()+
  #per replicate points
  geom_point(data=plotdata%>%
               filter(marker=="CO1",type%in%c("Autonomous_replicate","FAS_replicate")),
             aes(x=week,y=xval,col=type2),size=0.5, position=position_dodge(width=0.75),show.legend = FALSE)+
  
  #mean among replicates
  geom_errorbar(data=plotdata%>%
                  filter(marker=="CO1",type %in% c("Autonomous","FAS")),
                aes(x=week,y=xval,ymin=xval-se,ymax=xval+se,group=type2,col=type2),
                #linetype = "dotted",
                position=position_dodge(width=0.75),
                width=0.2,show.legend = FALSE)+
  geom_point(data=plotdata%>%filter(marker=="CO1",type %in% c("Autonomous","FAS")),
             position=position_dodge(width=0.75),size=3,
             aes(x=week,y=xval,group=type2,fill=type2),pch=21,col="black")+
  
  #shared species points and line
  geom_line(data=plotdata%>%
              filter(marker=="CO1",type=="Shared")%>%
              mutate(marker="CO1"),aes(x=week,y=xval,group=1),lty=2,col="grey80",lwd=0.6)+
  geom_point(data=plotdata%>%filter(marker=="CO1",type=="Shared"),
             aes(x=week,y=xval),shape=21,fill="grey80",size=1.5)+
  
  #total species 
  geom_line(data=plotdata%>%
              filter(marker=="CO1",type=="Total")%>%
              mutate(marker="CO1"),aes(x=week,y=xval,group=1),lty=2,col="grey30",lwd=0.6)+
  geom_point(data=plotdata%>%filter(marker=="CO1",type=="Total"),
             aes(x=week,y=xval),shape=21,fill="grey30",size=1.5)+
  
  # Labels for total and shared lines
  geom_text(data = plotdata %>%
              filter(marker == "CO1", type %in% c("Total", "Shared"), week == 9),
            aes(x = week, y = xval, label = type), 
            hjust = -0.1, vjust = -1.5, size = 3, color = c("grey80", "grey30"))  +
  
  #facet label
  facet_grid(marker~.)+
  
  #time of week labels
  geom_text(data=data.frame(week = c(1,9),
                            xval = rep(0,2), 
                            label = c("Aug 11 2023","Oct 16 2023")),aes(x=week,y=xval,label=label))+
  theme_bw()+
  labs(x="Sample week",y="Species detected",fill="")+
  theme(legend.position = "none",  
        strip.background = element_rect(fill="white"))+
  scale_y_continuous(limits = c(0,26))


#patchwork combination plot ----
p1 <- wrap_elements(p1_a/p1_b & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))) + 
  labs(tag = "Species detected") + # https://tidytales.ca/snippets/2022-12-22_patchwork-shared-axis-labels/
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )

ggsave("figures/marker_method_comparision.png",p1,width=8,height=5,units="in",dpi=300)


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
