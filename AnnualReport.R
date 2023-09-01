library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(nlme)
library(emmeans)
library(multcomp)
library(Hmisc)

options(contrasts=c("contr.helmert","contr.poly"))

#For the experiment, Audubon is only MM, EastKnoll and Lefty are only DM, and Trough is only WM

setwd("~/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/")

dat<-read.csv("Data/labdata2021.csv",stringsAsFactors=T)
head(dat)
dat$Community<-factor(dat$Community,levels=c("WM","SB","MM","DM","FF"))


###### Moisture - survey #####

dat2<-dat%>%
  filter(SoilMoistureP<150)%>%
  filter(Proj!="E")%>%
  filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one
  arrange(Site,Community)
dat2

ggplot(dat2, aes(x=Community,y=SoilMoistureP))+
  geom_boxplot()

m1<-dat2%>%
  group_by(Community)%>%
  summarise(mean=mean(SoilMoistureP), se=std.error(SoilMoistureP))
m1

pdf("Figs/MoistureSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Moisture - Experiment #####

dat2<-dat%>%
  filter(SoilMoistureP<150)%>%
  filter(Proj!="S")%>%
  arrange(Site,Community)%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_"),levels=c("DM_C","DM_E","MM_C","MM_E","WM_C","WM_E")))%>%
  mutate(Community=factor(Community,levels=c("DM","MM","WM")))
dat2

ggplot(dat2, aes(x=CommunityProj,y=SoilMoistureP))+
  geom_boxplot()

m1<-dat2%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(SoilMoistureP), se=std.error(SoilMoistureP))
m1

pdf("Figs/MoistureSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


###### Biomass - Survey #####

dat3<-dat%>%
  filter(is.na(Biomassg)==F)%>%
  filter(Proj!="E")%>%
  filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one
  arrange(Site,Community)

ggplot(dat3, aes(x=Community,y=Biomassg))+
  geom_boxplot()

m1<-dat3%>%
  group_by(Community)%>%
  summarise(mean=mean(Biomassg), se=std.error(Biomassg))
m1

pdf("Figs/BiomassSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Biomass (g)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Biomass - Experiment #####

dat3<-dat%>%
  filter(is.na(Biomassg)==F)%>%
  filter(Proj!="S")%>%
  arrange(Site,Community)%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_"),levels=c("DM_C","DM_E","MM_C","MM_E","WM_C","WM_E")))%>%
  mutate(Community=factor(Community,levels=c("WM","MM","DM")))

ggplot(dat3, aes(x=CommunityProj,y=Biomassg))+
  geom_boxplot()

m1<-dat3%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(Biomassg), se=std.error(Biomassg))
m1

pdf("Figs/BiomassSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Biomass (g)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+  
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()



###### AMF Survey ####

amf<-read.csv("Data/amf.csv",stringsAsFactors=T)
head(amf)
amf$Community<-factor(amf$Community,levels=c("WM","SB","MM","DM","FF"))

amf2<-amf%>%
  filter(Person!="KANELY")%>%
  filter(Proj!="E")%>%
  #filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one. i took this out b/c there were so few WM
  arrange(Site,Community)

ggplot(amf2, aes(x=Community,y=ArbusculesP))+
  geom_boxplot()
ggplot(amf2, aes(x=Community,y=VesiclesP))+
  geom_boxplot()
ggplot(amf2, aes(x=Community,y=HyphaeP))+
  geom_boxplot()

#Arbuscules
m1<-amf2%>%
  group_by(Community)%>%
  summarise(mean=mean(ArbusculesP), se=std.error(ArbusculesP))
m1

mod1<-lme(ArbusculesP~Community,random=~1|Site,data=amf2)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Community)))
anova(mod1, type="marginal")
summary(glht(mod1, linfct = mcp(Community = "Tukey")))

pdf("Figs/ArbusculesSurvey.pdf",width=2.2,height=2.2)
ggplot(data=mod1emm, aes(x=Community, y=emmean))+   
  geom_errorbar(aes(ymax=emmean+SE,ymin=emmean-SE),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Arbuscules (%)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()

#Vesicles
m1<-amf2%>%
  group_by(Community)%>%
  summarise(mean=mean(VesiclesP), se=std.error(VesiclesP))
m1

mod1<-lme(VesiclesP~Community,random=~1|Site,data=amf2)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Community)))
anova(mod1, type="marginal")
summary(glht(mod1, linfct = mcp(Community = "Tukey")))

pdf("Figs/VesiclesSurvey.pdf",width=2.2,height=2.2)
ggplot(data=mod1emm, aes(x=Community, y=emmean))+   
  geom_errorbar(aes(ymax=emmean+SE,ymin=emmean-SE),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Vesicles (%)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### AMF Experiment ####

#Arbuscules

amf3<-amf%>%
  filter(Person!="KANELY")%>%
  filter(Proj!="S")%>% #filter the survey plots
  filter(Site!="T")%>%
  arrange(Site)%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_")))
amf3$Proj<-recode_factor(amf3$Proj,"C"="Control","E"="Treatment","S"="Survey")

ggplot(amf3, aes(x=Proj,y=ArbusculesP))+
  geom_boxplot()
ggplot(amf3, aes(x=Proj,y=VesiclesP))+
  geom_boxplot()
ggplot(amf3, aes(x=Proj,y=HyphaeP))+
  geom_boxplot()

m1<-amf3%>%
  #group_by(Site,Community,Proj)%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(ArbusculesP), se=std.error(HyphaeP))
as.data.frame(m1)

pdf("Figs/ArbusculesExperiment.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean,group=Community))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Arbuscules (%)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


mod1<-lme(ArbusculesP~Proj+Community+Proj*Community,random=~1|Site/Plot,data=amf3)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Proj|Community)))
anova(mod1, type="marginal")

mod1<-lme(ArbusculesP~CommunityProj,random=~1|Site/Plot,data=amf3)
summary(glht(mod1, linfct = mcp(CommunityProj = "Tukey")))
summary(glht(mod1, linfct = mcp(CommunityProj=c("DM_E-DM_C=0","MM_E-MM_C=0"))))

pdf("Figs/AMF.pdf",width=2.2,height=2.2)
ggplot(data=mod1emm, aes(x=Proj, y=emmean))+   
  geom_errorbar(aes(ymax=emmean+SE,ymin=emmean-SE),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Arbuscules (%)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()

