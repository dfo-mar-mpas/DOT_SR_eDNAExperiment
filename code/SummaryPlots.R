## Plotting DOT comparison data -----

#load libraries ---
library(tidyverse)
library(patchwork)
library(vegan)
library(ggnewscale)
library(vegan)
library(drc)

#functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Function to estimate number of sites needed to achieve estimated richness
estimate_sites <- function(fit, richness_estimate) {
  coef <- coef(fit)
  asymptote <- coef[1]
  rate <- coef[2]
  return((richness_estimate * rate) / (asymptote - richness_estimate))
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

#summary stats for paper
mod <- lm(nspecies ~ week*type,data=dot_df_rich_rep%>%filter(marker=="12S")) #'week' as a surrogate of date. No significant trends but can cross refernce if needed. 
anova(mod)

mod <- lm(nspecies ~ week*type,data=dot_df_rich_rep%>%filter(marker=="CO1"))
anova(mod)

mod <- aov(nspecies ~ type,data=dot_df_rich_rep%>%filter(marker=="12S"))
anova(mod)

mod <- aov(nspecies ~ type,data=dot_df_rich_rep%>%filter(marker=="CO1"))
anova(mod)

#how much species per marker and method
dot_df_rich_rep%>%
  group_by(marker,type)%>%
  summarize(mn=mean(nspecies),
            sd=sd(nspecies))%>%
  ungroup()%>%
  data.frame()

dot_df_rich_rep%>%
  group_by(marker)%>%
  summarize(mn=mean(nspecies),
            sd=sd(nspecies))%>%
  ungroup()%>%
  data.frame()

#how much overlap in species composition was there per marker
dot_df_week_rich_overlap%>%
  group_by(marker)%>%
  summarize(mn=mean(overlap),
            sd=sd(overlap))%>%
  ungroup()%>%
  data.frame()

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
             aes(x=week,y=xval),shape=21,fill="grey30",size=1.5)+
  
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


#eDNA Concentrations ------
conc_data <- read.csv("data/eDNA_Concentrations.csv")%>%
             mutate(type=ifelse(Method == "DOT","Autonomous","FAS"),
                    concentration=as.numeric(ifelse(Concentration == "<0.05",0.05,Concentration)),
                    rep=case_when(DOT_unit==1003 ~ 1,
                                  DOT_unit==1004 ~ 2,
                                  DOT_unit==1005 ~ 3,
                                  TRUE ~ NA),
                    rep=as.character(rep))


#compare the difference in eDNA concentration as a function of time

comp_df <- conc_data%>%
           filter(DurationPreserved %in% intersect(conc_data%>%filter(type=="Autonomous")%>%pull(DurationPreserved)%>%unique(),
                                                   conc_data%>%filter(type=="FAS")%>%pull(DurationPreserved%>%unique())))%>%
           rename(duration=DurationPreserved)%>%
           select(rep,type,duration,concentration)

ggplot(comp_df,aes(x=duration,y=concentration,col=interaction(rep,type),group=interaction(rep,type)))+
  #facet_wrap(~type,nrow=2)+
  geom_point()+
  geom_line()+
  theme_bw()+
  stat_smooth(aes(group = type,col=type),method="lm")

#standardsize the data for the preservation comparison
comp_df_stand <- comp_df%>%
                 group_by(type,rep)%>%
                 mutate(stand = concentration/max(concentration))%>%
                 ungroup()%>%
                 data.frame()

# Fit the linear model
model <- lm(stand ~ duration * type, data = data)

# Summarize the model
summary(model)

# Generate ANOVA table
anova_table <- anova(model)

#Set up prediction data for the plot
pred_data <- data %>%
  group_by(type) %>%
  reframe(duration = seq(min(duration), max(duration), length.out = 100))

pred_data <- pred_data %>%
  mutate(pred = predict(model, newdata = pred_data, interval = "confidence"))

pred_data$stand = pred_data[,3]

plotcols <- c("cornflowerblue","tomato3")
  
  p1 <- ggplot(data, aes(x = duration, color = type)) +
    geom_point(aes(y=stand)) +
    geom_line(data = pred_data, aes(x = duration, y = pred[,1], color = type), size = 1) +
    geom_ribbon(data = pred_data, aes(x = duration, y = pred[,1], ymin = pred[,2], ymax = pred[,3],fill=type), alpha = 0.2) +
    geom_smooth(aes(y=stand,group = interaction(type, rep)), method = "lm", se = FALSE, linetype = "dashed", size = 0.5) +
    theme_bw() +
    scale_fill_manual(values=plotcols)+
    scale_colour_manual(values=plotcols)+
    labs(x="Preservation duration (days)",y="Standardized DNA yeild",col="",fill="")+
    theme(legend.position = c(0, 0),  # Position legend in the top left
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.title = element_blank(),
          strip.background = element_rect(fill="white"))
  
  ggsave("figures/DNAYield_time.png",p1,height=6,width=8,units="in",dpi=300)


  #alternative test - does the differnce between the two methods vary as a function time. 
  
comp_df2 <- comp_df%>%
            spread(type,concentration)%>%
            mutate(diff=Autonomous - FAS)

model2 <- lm(diff~duration,data=comp_df2)
anova(model2) #not a significant slope. 

ggplot(data=comp_df2,aes(x=duration,y=diff))+
         geom_point(aes(col=rep,group=rep))+
         geom_line(aes(col=rep,group=rep))+
         geom_hline(yintercept=0,lty=2)+
         stat_smooth(method="lm")+
         theme_bw()+
         theme(legend.position = "none")


### NMDS plots ------------------

#data setup

#12S marker
community_matrix_12S <- dot_data_raw%>%
                        select(-2)%>%
                        gather("sample","reads",2:48)%>%
                        pivot_wider(names_from = BLAST,values_from = reads)%>%
                        mutate(sample2 = sample,
                               sample = gsub("X","Autonomous_",sample))%>%
                        as_tibble()%>%
                        separate(sample2,c("week","name"))%>%
                        mutate(sample3 = gsub("Parallel","",name),
                               sample4 = case_when(grepl("Parallel",sample) ~ paste0("FAS_",sample3,"_",week),
                                                   grepl("Control",sample) ~ paste0("FAS_",1001,"_",week),
                                                   TRUE~sample))%>%
                        select(-c(sample3,sample,week,name))%>%
                        rename(sample=sample4)%>%
                        column_to_rownames(var="sample")%>%
                        mutate(tot = rowSums(.,na.rm=T))%>%
                        filter(tot>0)%>% #filter out any zero catches
                        select(-tot)%>%
                        decostand(method = "total")

groups_12S <- c(rep("Autonomous",sum(grepl("Autonomous",rownames(community_matrix_12S)))),
                rep("FAS",sum(!grepl("Autonomous",rownames(community_matrix_12S)))))


#COI Marker
community_matrix_CO1 <- dot_data_raw_invert%>%
                        select(-c(1:4,6))%>%
                        mutate(Species = gsub("_"," ",Species))%>%
                        gather("sample","reads",2:48)%>%
                        pivot_wider(names_from = Species,values_from = reads)%>%
                        mutate(sample_org = sample,
                               sample = gsub("DOT","Autonomous_",sample))%>%
                        as_tibble()%>%
                        separate(sample_org,c("type","week"))%>%
                        mutate(sample=case_when(grepl("SR1",type) ~ paste0("FAS_",1001,"_","Wk",week),
                                                grepl("SR3",type) ~ paste0("FAS_",1003,"_","Wk",week),
                                                grepl("SR4",type) ~ paste0("FAS_",1004,"_","Wk",week),
                                                grepl("SR5",type) ~ paste0("FAS_",1005,"_","Wk",week),
                                                TRUE ~ sample))%>%
                        select(-c(type,week))%>%
                        column_to_rownames(var="sample")%>%
                        mutate(tot = rowSums(.,na.rm=T))%>%
                        filter(tot>0)%>% #filter out any zero catches
                        select(-tot)%>%
                        decostand(method = "total")

groups_CO1 <- c(rep("Autonomous",sum(grepl("Autonomous",rownames(community_matrix_CO1)))),
                rep("FAS",sum(!grepl("Autonomous",rownames(community_matrix_CO1)))))

#Run the metaMDS analyses based on Bray Curtis and Euclidean distance

Bray_12s <-metaMDS(community_matrix_12S, distance="bray", k=12, trymax = 300, maxit=500)
Euclid_12s <- metaMDS(community_matrix_12S%>%dist(., method="binary"),distance="binary",k=12, trymax = 300, maxit=500)

Bray_CO1 <-metaMDS(community_matrix_CO1, distance="bray", k=12, trymax = 300, maxit=500)
Euclid_CO1 <- metaMDS(community_matrix_CO1%>%dist(., method="binary"),distance="binary",k=12, trymax = 300, maxit=500)

#data processing -------------
data.scores <- as.data.frame(scores(Bray_12s,"sites"))%>%
              mutate(id=rownames(.),
                     method="Bray-Curtis",
                     stress=round(Bray_12s$stress,3))%>%
              separate(id,c("type","rep","week"),sep="_")%>%
              rbind(.,
                    as.data.frame(scores(Euclid_12s,"sites"))%>%
                      mutate(id=rownames(.),
                             method="Euclidean distance",
                             stress=round(Euclid_12s$stress,3))%>%
                      separate(id,c("type","rep","week"),sep="_"))%>%
              select(NMDS1,NMDS2,method,type,rep,week,stress)%>%
              mutate(marker="12S")%>%
              rbind(.,as.data.frame(scores(Bray_CO1,"sites"))%>%
                      mutate(id=rownames(.),
                             method="Bray-Curtis",
                             stress=round(Bray_CO1$stress,3))%>%
                      separate(id,c("type","rep","week"),sep="_")%>%
                      rbind(.,
                              as.data.frame(scores(Euclid_CO1,"sites"))%>%
                              mutate(id=rownames(.),
                                     method="Euclidean distance",
                                     stress=round(Euclid_CO1$stress,3))%>%
                              separate(id,c("type","rep","week"),sep="_"))%>%
                              select(NMDS1,NMDS2,method,type,rep,week,stress)%>%
                              mutate(marker="CO1"))
                               
species.scores <- as.data.frame(scores(Bray_12s,"species"))%>%
                  mutate(method = "Bray-Curtis")%>%
                  rbind(.,as.data.frame(scores(Euclid_12s,"species"))%>%
                          mutate(method="Euclidean distance"))%>%
                  mutate(marker="12S")%>%rbind(.,
                                               as.data.frame(scores(Bray_CO1,"species"))%>%
                                                 mutate(method = "Bray-Curtis")%>%
                                                 rbind(.,as.data.frame(scores(Euclid_CO1,"species"))%>%
                                                         mutate(method="Euclidean distance"))%>%
                                                 mutate(marker="CO1"))%>%
                  rownames_to_column(var="species")%>%
                  select(species,NMDS1,NMDS2,marker)
                                                  

hull_data <- data.scores%>%
            group_by(marker,method,type)%>%
            slice(chull(x=NMDS1,y=NMDS2))%>%
            mutate(weekn = as.numeric(gsub("Wk","",week)))


## make plots

nmds_plot <- ggplot(data=hull_data)+
              #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),size=2,col="grey70")+
              geom_polygon(aes(x=NMDS1,y=NMDS2,fill=type),alpha=0.30)+
              labs(fill="")+
              geom_point(aes(x=NMDS1,y=NMDS2,fill=type),col="black",shape=21,size=1)+
              geom_point(aes(x=NMDS1,y=NMDS2,fill=type,size=log10(weekn+1)),col="black",shape=21,show.legend = FALSE)+
              scale_size_continuous(range = c(2, 7),breaks=log10(c(1,2,5,9)+1),labels=c(1,2,5,9))+
              theme_bw()+
              facet_grid(method~marker,scales="free")+
              theme(strip.background = element_rect(fill="white"),
                    legend.position="inside",
                    legend.position.inside = c(0.085,0.96),
                    legend.background = element_blank())
              
                            
                                
ggsave("figures/nmds_plot_all.png",nmds_plot,height=6,width=8,units="in",dpi=300)


# peform a multivariate analysis of variance. 

## 12s
# Extract sample identifiers
samples <- rownames(community_matrix_12S)

# Extract treatments and weeks from sample identifiers
treatment <- ifelse(grepl("Autonomous", samples), "Autonomous", "FAS")
week <- as.numeric(sub(".*Wk([0-9]+).*", "\\1", samples))
replicate_unit <- samples%>%
                  data.frame()%>%
                  rename(sample = 1)%>%
                  separate(sample,c("treament","replicate","week"))%>%
                  pull(replicate)

# Create an environmental data frame
env_data <- data.frame(treatment = treatment, week = week,replicate_unit = replicate_unit)

# Conduct PERMANOVA
adonis_12s_bray <- adonis2(community_matrix_12S ~ treatment * week + replicate_unit, data = env_data, method = "bray")
adonis_12s_euclid <- adonis2(community_matrix_12S ~ treatment * week, data = env_data, method = "euclidean")

# Display the results
print(adonis_12s_bray)
print(adonis_12s_euclid)

## COI
# Extract sample identifiers
samples_CO1 <- rownames(community_matrix_CO1)

# Extract treatments and weeks from sample identifiers
treatment_CO1 <- ifelse(grepl("Autonomous", samples_CO1), "Autonomous", "FAS")
week_CO1 <- as.numeric(sub(".*Wk([0-9]+).*", "\\1", samples_CO1))
replicate_unit_CO1 <- samples_CO1%>%
                  data.frame()%>%
                  rename(sample = 1)%>%
                  separate(sample,c("treament","replicate","week"))%>%
                  pull(replicate)


# Create an environmental data frame
env_data_CO1 <- data.frame(treatment = treatment_CO1, week = week_CO1,replicate_unit = replicate_unit_CO1)

# Conduct PERMANOVA
adonis_CO1_bray <- adonis2(community_matrix_CO1 ~ treatment * week + replicate_unit, data = env_data, method = "bray")
adonis_CO1_euclid <- adonis2(community_matrix_CO1 ~ treatment * week + replicate_unit, data = env_data, method = "euclidean")

# Display the results
print(adonis_CO1_bray)
print(adonis_CO1_euclid)


# Species accumulation curves ------

meta_12s <- rownames(community_matrix_12S)%>%
            data.frame()%>%
            rename(sample = 1)%>%
            separate(sample,c("treament","replicate","week"))

meta_CO1 <- rownames(community_matrix_CO1)%>%
            data.frame()%>%
            rename(sample = 1)%>%
            separate(sample,c("treament","replicate","week"))

 
SA_df_12S_auto <- community_matrix_12S[meta_12s == "Autonomous",]%>%
                  select_if(~ sum(.) != 0)

SA_df_12S_FAS <- community_matrix_12S[meta_12s == "FAS",]%>%
                  select_if(~ sum(.) != 0)

SA_df_CO1_auto <- community_matrix_CO1[meta_CO1 == "Autonomous",]%>%
                  select_if(~ sum(.) != 0)

SA_df_CO1_FAS <- community_matrix_CO1[meta_CO1 == "FAS",]%>%
                select_if(~ sum(.) != 0)

#species accumulation curves
sa_autonomous_12S <- specaccum(SA_df_12S_auto, method = "random")
sa_fas_12S <- specaccum(SA_df_12S_FAS, method = "random")

sa_autonomous_CO1 <- specaccum(SA_df_CO1_auto, method = "random")
sa_fas_CO1 <- specaccum(SA_df_CO1_FAS, method = "random")

#estimate richness
richness_est_12S <- specpool(SA_df_12S_auto)%>%
                       data.frame()%>%
                       mutate(method="Autonomous")%>%
                       rbind(.,specpool(SA_df_12S_FAS)%>%
                               data.frame()%>%
                               mutate(method="FAS"))%>%
                       mutate(marker="12S")

richness_est_CO1 <- specpool(SA_df_CO1_auto)%>%
                  data.frame()%>%
                  mutate(method="Autonomous")%>%
                  rbind(.,specpool(SA_df_CO1_FAS)%>%
                          data.frame()%>%
                          mutate(method="FAS"))%>%
                  mutate(marker="CO1")


## Fit Michaelis-Menten model to the accumulation curves
fit_autonomous_12S <- drm(sa_autonomous_12S$richness ~ sa_autonomous_12S$sites, fct = MM.2())
fit_fas_12S <- drm(sa_fas_12S$richness ~ sa_fas_12S$sites, fct = MM.2())

fit_autonomous_CO1 <- drm(sa_autonomous_CO1$richness ~ sa_autonomous_CO1$sites, fct = MM.2())
fit_fas_CO1 <- drm(sa_fas_CO1$richness ~ sa_fas_CO1$sites, fct = MM.2())

#dataframe for the estimated assymptotic diversity
sites_needed_autonomous_12S <- estimate_sites(fit_autonomous_12S, richness_est_12S%>%filter(method=="Autonomous")%>%pull(chao))
sites_needed_fas_12S <- estimate_sites(fit_fas_12S, richness_est_12S%>%filter(method=="FAS")%>%pull(chao))

sites_needed_autonomous_CO1 <- estimate_sites(fit_autonomous_CO1, richness_est_CO1%>%filter(method=="Autonomous")%>%pull(chao))
sites_needed_fas_CO1 <- estimate_sites(fit_fas_CO1, richness_est_CO1%>%filter(method=="FAS")%>%pull(chao))

#add to the richness dfs
richness_est_12S$sites <- c(sites_needed_autonomous_12S,sites_needed_fas_12S)
richness_est_12S$msite <- sa_12S_df%>%group_by(method)%>%summarize(max=max(Sites))%>%ungroup()%>%pull(max)
richness_est_CO1$sites <- c(sites_needed_autonomous_CO1,sites_needed_fas_CO1)
richness_est_CO1$msite <- sa_CO1_df%>%group_by(method)%>%summarize(max=max(Sites))%>%ungroup()%>%pull(max)

# Convert to data frames for ggplot2
sa_12S_df <- data.frame(Sites = sa_autonomous_12S$sites, 
                            Richness = sa_autonomous_12S$richness, 
                            SD = sa_autonomous_12S$sd)%>%
                    mutate(method="Autonomous")%>%
  rbind(., data.frame(Sites = sa_fas_12S$sites, 
                     Richness = sa_fas_12S$richness, 
                     SD = sa_fas_12S$sd)%>%
             mutate(method="FAS"))%>%
  mutate(marker="12S")

sa_CO1_df <- data.frame(Sites = sa_autonomous_CO1$sites, 
                        Richness = sa_autonomous_CO1$richness, 
                        SD = sa_autonomous_CO1$sd)%>%
  mutate(method="Autonomous")%>%
  rbind(., data.frame(Sites = sa_fas_CO1$sites, 
                      Richness = sa_fas_CO1$richness, 
                      SD = sa_fas_CO1$sd)%>%
          mutate(method="FAS"))%>%
  mutate(marker="CO1")

plot_df <- rbind(sa_12S_df,sa_CO1_df)
plot_df_rich <- rbind(richness_est_12S,richness_est_CO1)


sa_plot <- ggplot(data=plot_df,aes(group=method,fill=method)) +
          geom_line(data=plot_df,aes(color=method,x = Sites, y = Richness),lwd=1.5) +
          geom_ribbon(data=plot_df,aes(color=method,x = Sites, ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2)+
          geom_errorbar(data=plot_df_rich,aes(x=rep(c(28,29),2),ymin = chao-chao.se,ymax=chao+chao.se),size=0.5)+
          geom_point(data=plot_df_rich,aes(x=rep(c(28,29),2),y = chao),pch=21,size=2)+
          theme_bw()+
          facet_wrap(~marker,ncol=2)+
          labs(x = "Number of Samples", y = "Species Richness",fill="",color="")+
          theme(strip.background = element_rect(fill="white"),
                legend.position="inside",
                legend.position.inside = c(0.085,0.93),
                legend.title=element_blank(),
                legend.background = element_blank())

ggsave("figures/species_accum_plots.png",sa_plot,width=8,height=5,units="in",dpi=300)               
    

+
  labs(title = "Species Accumulation Curves", x = "Number of Samples", y = "Species Richness") +
  theme_minimal() +
  scale_color_manual(values = c("Autonomous" = "blue", "FAS" = "red")) +
  guides(fill = guide_legend(title = "Sampling Method"))

