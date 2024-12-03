#LTER data analysis


##### Chiara's 2020 itex data #####

moist<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/LTER/black_sand_soil_moisture_2020.cf.data.csv",stringsAsFactors = T)
head(moist)
moist$local_site<-factor(moist$local_site,levels=c("Trough","Audubon","Lefty","East_Knoll"))
moist$treat<-paste(moist$snow_treatment,moist$warm_treatment,sep="")

moist2<-moist%>%
  filter((snow_treatment=="Control"&warm_treatment=="Control")|(snow_treatment=="Early"&warm_treatment=="OTC"))%>%
  group_by(local_site,treat,snow_treatment,warm_treatment,date)%>%
  summarise(mean=mean(VWC),se=std.error(VWC))
data.frame(moist2)

ggplot(data=moist2, aes(x=treat, y=mean,color=date))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  theme_classic()+
  #theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  facet_wrap(~local_site,nrow=1)


moist2<-moist%>%
  filter((snow_treatment=="Control"&warm_treatment=="Control")|(snow_treatment=="Early"&warm_treatment=="OTC"))%>%
  group_by(local_site,treat,snow_treatment,warm_treatment)%>%
  summarise(mean=mean(VWC),se=std.error(VWC))

ggplot(data=moist2, aes(x=treat, y=mean,color=local_site))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  theme_classic()


moist2<-moist%>%
  group_by(local_site,treat,snow_treatment,warm_treatment)%>%
  summarise(mean=mean(VWC),se=std.error(VWC))

ggplot(data=moist2, aes(x=treat, y=mean,color=local_site))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  theme_classic()



##### 2018 only control and black sand data #####

moistbs<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/LTER/black_sand_soil_moisture_2018.cf.data.csv",stringsAsFactors = T)
head(moistbs)
moistbs$local_site<-factor(moistbs$local_site,levels=c("Trough","Audubon","Lefty","East Knoll","Soddie"))

moistbs2<-moistbs%>%
  filter(local_site!="Soddie")%>%
  group_by(local_site,snow_treatment,time_period)%>%
  summarise(mean=mean(VWC),se=std.error(VWC))

ggplot(data=moistbs2, aes(x=snow_treatment, y=mean,color=local_site))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  theme_classic()+
  facet_wrap(~time_period)

ggplot(data=moistbs2, aes(x=snow_treatment, y=mean,color=time_period))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  theme_classic()+
  facet_wrap(~local_site,nrow=1)


##### Veg data #####
library(cluster)
library(dendextend)
#for clustering see https://www.datacamp.com/tutorial/hierarchical-clustering-R

veg<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/LTER/gsle_itex_subplots.cf.data.csv",stringsAsFactors = T)
veg$treat<-paste(veg$snow_treatment,veg$warm_treatment,sep="")
veg$SiteTreatPlot<-paste(veg$local_site,veg$snow_treatment,veg$subplot,sep="")
veg$SitePlot<-paste(veg$local_site,veg$subplot,sep="")
veg$species2<-paste(veg$NWT_code,veg$former_ID,sep="")#merging nwt_code adn former_id so that carex specie can be distinguished, it doesn't change anything

veg

veg2<-veg%>%
  filter(year==2018)%>%
  filter(treat%in%c("ControlControl","EarlyOTC"))%>%
  dplyr::select(SiteTreatPlot,local_site,species2,hits)%>% #NWT_code
  pivot_wider(names_from = species2,values_from = hits,values_fn=sum,values_fill = 0) #names_from = NWT_code
  #dplyr::summarise(n = dplyr::n(), .by = c(local_site, subplot,year, treat, NWT_code))
env2<-veg2
veg2<-data.frame(veg2)
rownames(veg2)<-veg2$SiteTreatPlot
veg2<-veg2[,-c(1:2)]
rowSums(veg2)

dist_mat <- vegdist(veg2, method = 'bray')
clus<-agnes(dist_mat,method = 'flexible',par.method = 0.625)# par.method = 0.625
clus.h <- as.hclust(clus)
plot(clus.h)

clus.h$order
groups <- cutree(clus.h, k = 2)

avg_dend_obj <- as.dendrogram(clus.h)
avg_col_dend <- color_branches(avg_dend_obj, k = 2)
plot(avg_col_dend)


vegmds<-metaMDS(veg2,distance="bray")#autotransform=F
vegmdsscores<-scores(vegmds)$sites
env2<-cbind(env2,vegmdsscores)

vegmds<-capscale(veg2~env2$local_site,distance="bray")#autotransform=F
vegmdsscores<-scores(vegmds)$sites
#cbind(as.character(env2$local_site),rownames(vegmdsscores))
env2<-cbind(env2,vegmdsscores)

vegmds<-capscale(veg2~1,distance="bray")#autotransform=F
vegmdsscores<-scores(vegmds)$sites
#cbind(as.character(env2$local_site),rownames(vegmdsscores))
env2<-cbind(env2,vegmdsscores)
env2

ggplot(data=env2, aes(x=CAP1, y=CAP2,color=local_site))+   
  geom_point(size=1.8,show.legend = T)+
  theme_classic()
ggplot(data=env2, aes(x=NMDS1, y=NMDS2,color=local_site))+   
  geom_point(size=1.8,show.legend = T)+
  theme_classic()
ggplot(data=env2, aes(x=MDS1, y=MDS2,color=local_site))+   
  geom_point(size=1.8,show.legend = T)+
  theme_classic()
