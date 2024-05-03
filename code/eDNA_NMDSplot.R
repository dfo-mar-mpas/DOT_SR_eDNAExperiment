library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)
library(vroom)

specs <- read.csv(file = "data/12S_Taxonomy_Table_FilteredFinal_grouped.csv", header = T)

#get unique species
sort(unique(specs$Sp))
#42 species

commat <- specs[,3:48]
groups<- c(rep("DOT",27), rep("SmithRoot",19))

#Transpose the table for vegan
commat2<-t(commat)
colnames(commat2)<-specs[,1] #this selects the 3rd column from our taxon table and inserts the species IDs into this matrix

#commat3<-commat2 %>% group_by(colnames())


dot.nmds <-metaMDS(commat2, distance="bray", k=12, trymax = 200, maxit=500)
plot(dot.nmds) #this is not very informative without labels!

#slightly better plot here
ordiplot(dot.nmds, type='n')
orditorp(dot.nmds, display="species", col="red", air=0.01)
orditorp(dot.nmds, display="sites", cex=1, air=0.01)



#plot NMDS with ggplot to look nicer
data.scores <- as.data.frame(scores(dot.nmds,"sites")) %>%
  mutate(ID=rownames(.), group=groups)
species.scores <- as.data.frame(scores(dot.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(group) %>%
  slice(chull(x=NMDS1,y=NMDS2))

nmdstheme <- theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))


# regular ggplot
p1 <- ggplot() +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,color=group, fill=group),alpha=0.20) + # add the hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels or remove this if too messy
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=group),size=3) +
  #scale_colour_manual()
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2, label=ID), size=6)+
  #scale_color_gradient(low="blue", high="black")+
  scale_colour_brewer(palette = "Accent") +
  coord_equal()+
  #scale_shape_manual(values=c(15,7,18,16,10,8,17))+
  geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(dot.nmds$stress,3),"k =",dot.nmds$ndim)))+
  nmdstheme
p1


ggsave(filename = "DOT_NMDS.png", plot = p1, device = "png", path = "figures/", width = 16, height=12, units = "in", dpi=400, bg = "white")


#Species accumulation curve and other info 
gg<- specaccum(commat2, method="exact", ci.type="polygon")
ggg <-data.frame(gg$sites,gg$richness,gg$sd)

ggplot()+geom_line(data=ggg, aes(x=gg.sites,y=gg.richness),colour="navyblue", lwd=1.25)+
  geom_errorbar(data=ggg,aes(x=gg.sites,ymin=gg.richness-gg.sd, ymax=gg.richness+gg.sd), width=0, colour="navyblue")+
  ylab(label = "Species richness")+
  xlab(label = "# Samples")+
  theme_bw()+
  theme(text=element_text(size=20))
ggsave(filename = "DOT_specaccum.png", plot = last_plot(), device = "png", path = "figures/", width = 14, height=10, units = "in", dpi=400, bg = "white")

#rarefaction
trask.rare<-rarecurve(commat2, step=10)


#this gives us Shannon index values per sample 
shan <- data.frame(diversity(commat2, index="shannon"))
#we can aggregate these so we get a site level index
mean.shan <- aggregate(. ~ substr(rownames(shan), 1, 3), shan, mean, na.rm = TRUE)


########################################################################################
##Now the COI data
##################

coi <-read.csv("data/DOTCOI_feature_table_export_filtered.csv")
#36 x 53 dataframe 
coi.commat <- coi[,7:length(colnames(coi))]
groups<- c(rep("DOT",27), rep("SmithRoot",20))


#Transpose the table for vegan
coi.commat2<-t(coi.commat)
colnames(coi.commat2)<-coi$Species 



dot.nmds.coi <-metaMDS(coi.commat2, distance="bray", k=12, trymax = 200, maxit=500)
plot(dot.nmds.coi) #this is not very informative without labels!

#slightly better plot here
ordiplot(dot.nmds.coi, type='n')
orditorp(dot.nmds.coi, display="species", col="red", air=0.01)
orditorp(dot.nmds.coi, display="sites", cex=1, air=0.01)



#plot NMDS with ggplot to look nicer
data.scores.coi <- as.data.frame(scores(dot.nmds.coi,"sites")) %>%
  mutate(ID=rownames(.), group=groups)
species.scores.coi <- as.data.frame(scores(dot.nmds.coi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores.coi$species <- rownames(species.scores.coi) 

hull.data.coi <- data.scores.coi %>%
  as.data.frame() %>%
  group_by(group) %>%
  slice(chull(x=NMDS1,y=NMDS2))


# regular ggplot
p2 <- ggplot() +
  geom_polygon(data=hull.data.coi,aes(x=NMDS1,y=NMDS2,color=group, fill=group),alpha=0.20) + # add the hulls
  geom_text(data=species.scores.coi,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels or remove this if too messy
  geom_point(data=data.scores.coi,aes(x=NMDS1,y=NMDS2,colour=group),size=3) +
  #scale_colour_manual()
  geom_text(data=data.scores.coi,aes(x=NMDS1,y=NMDS2, label=ID), size=6)+
  #scale_color_gradient(low="blue", high="black")+
  scale_colour_brewer(palette = "Accent") +
  coord_equal()+
  #scale_shape_manual(values=c(15,7,18,16,10,8,17))+
  geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(dot.nmds.coi$stress,3),"k =",dot.nmds.coi$ndim)))+
  nmdstheme
p2


ggsave(filename = "DOT_NMDS_COI.png", plot = p2, device = "png", path = "figures/", width = 16, height=12, units = "in", dpi=400, bg = "white")


#Species accumulation curve and other info 
gg<- specaccum(coi.commat2, method="exact", ci.type="polygon")
ggg <-data.frame(gg$sites,gg$richness,gg$sd)

ggplot()+geom_line(data=ggg, aes(x=gg.sites,y=gg.richness),colour="navyblue", lwd=1.25)+
  geom_errorbar(data=ggg,aes(x=gg.sites,ymin=gg.richness-gg.sd, ymax=gg.richness+gg.sd), width=0, colour="navyblue")+
  ylab(label = "Species richness")+
  xlab(label = "# Samples")+
  theme_bw()+
  theme(text=element_text(size=20))
ggsave(filename = "DOT_specaccum_COI.png", plot = last_plot(), device = "png", path = "figures/", width = 14, height=10, units = "in", dpi=400, bg = "white")

#rarefaction
coi.rare<-rarecurve(coi.commat2, step=10)


#this gives us Shannon index values per sample 
shan <- data.frame(diversity(commat2, index="shannon"))
#we can aggregate these so we get a site level index
mean.shan <- aggregate(. ~ substr(rownames(shan), 1, 3), shan, mean, na.rm = TRUE)














save.image("data/DOTeDNA_CommunityStats.RData")

