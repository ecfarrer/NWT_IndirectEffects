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


#### 2023 Final report ####


labdat<-read.csv("Data/Niwot_IndirectEffects_2021_LabData.csv",stringsAsFactors = T)
labdat2<-labdat%>%
  filter(SoilMoisturePercent>0)%>%
  filter(!(SampleID=="ES-FF2"&Replicate==1))%>%
  filter(!(SampleID=="AS-DM3"&Replicate==1))%>%
  filter(!(SampleID=="LE-DM1c"&Replicate==1))%>%
  filter(!(SampleID=="AE-MM1a"&Replicate==1))%>%
  filter(!(SampleID=="TC-WM1a"&Replicate==1))%>%
  filter(!(SampleID=="TE-WM3c"&Replicate==1))%>%
  dplyr::select(SampleID:Replicate,SoilMoisturePercent,Biomassg)%>%
  group_by(SampleID,SiteProj,Site,Proj,Community,Plot,PlantID)%>%
  summarise(across(SoilMoisturePercent:Biomassg,~mean(.x,na.rm=T)))%>%
  arrange(SampleID)
labdat2<-data.frame(labdat2)
labdat2$Community<-factor(labdat2$Community,levels=c("SB","WM","MM","DM","FF"))
View(labdat2)

fielddat<-read.csv("Data/Niwot_IndirectEffects_2021_FieldData.csv",stringsAsFactors = T)
fielddat$Community<-factor(fielddat$Community,levels=c("SB","WM","MM","DM","FF"))

dim(labdat2)
dim(fielddat)

##### Getting survey data together #####
labdat2S<-labdat2%>%
  filter(Proj!="E")%>%
  #filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one
  arrange(Site,Community)
labdat2S

fielddatS<-fielddat%>%
  filter(Proj!="E")%>%
  arrange(Site,Community)

##### Getting Experiment data together #####
labdat2E<-labdat2%>%
  filter(Proj!="S")%>%
  arrange(Site,Community)%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_"),levels=c("WM_C","WM_E","MM_C","MM_E","DM_C","DM_E")))%>%
  mutate(Community=factor(Community,levels=c("WM","MM","DM")))
labdat2E

fielddatE<-fielddat%>%
  filter(Proj!="S")%>%
  arrange(Site,Community)%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_"),levels=c("WM_C","WM_E","MM_C","MM_E","DM_C","DM_E")))%>%
  mutate(Community=factor(Community,levels=c("WM","MM","DM")))



###### Moisture - survey #####

ggplot(labdat2S, aes(x=Community,y=SoilMoisturePercent))+
  geom_boxplot()

m1<-labdat2S%>%
  group_by(Community)%>%
  summarise(mean=mean(SoilMoisturePercent), se=std.error(SoilMoisturePercent))
m1

pdf("Figs/MoistureSurveynew.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Moisture - Experiment #####

ggplot(labdat2E, aes(x=CommunityProj,y=SoilMoisturePercent))+
  geom_boxplot()

m1<-labdat2E%>%
  filter(!is.na(SoilMoisturePercent))%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(SoilMoisturePercent), se=std.error(SoilMoisturePercent))
m1
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/MoistureExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()



###### Biomass - survey #####
#labdat2Sb<-labdat2S%>%
#  filter(!is.na(Biomassg))

ggplot(labdat2S, aes(x=Community,y=Biomassg))+
  geom_boxplot()

labdat2S%>%filter(Biomassg>2)

#for now I'm removing the three huge plants that had biomass >2g. These were huge, they had 7, 10 and 10 flowering stalks, all in fell field. ES-FF4, ES-FF6, TS-FF6

labdat2Sb<-labdat2S%>%
  filter(Biomassg<2)

m1<-labdat2Sb%>%
  group_by(Community)%>%
  summarise(mean=mean(Biomassg,na.rm=T), se=std.error(Biomassg,na.rm=T))
m1

pdf("Figs/BiomassSurveynew.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Biomass (g)")+
  ylim(.4,1.5)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Biomass - Experiment #####

labdat2E

ggplot(labdat2E, aes(x=CommunityProj,y=Biomassg))+
  geom_boxplot()

m1<-labdat2E%>%
  filter(!is.na(Biomassg))%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(Biomassg), se=std.error(Biomassg))
m1
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/BiomassExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Biomass (g)")+
  theme_classic()+
  ylim(.4,1.5)+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


###### Leaf number - survey #####
fielddatS
ggplot(fielddatS, aes(x=Community,y=LeafNumber))+
  geom_boxplot()

m1<-fielddatS%>%
  group_by(Community)%>%
  summarise(mean=mean(LeafNumber,na.rm=T), se=std.error(LeafNumber,na.rm=T))
m1

pdf("Figs/LeafNumberSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylim(7,11)+
  ylab("Number of leaves")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Leaf number - Experiment #####

ggplot(fielddatE, aes(x=CommunityProj,y=LeafNumber))+
  geom_boxplot()

m1<-fielddatE%>%
  filter(!is.na(LeafNumber))%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(LeafNumber), se=std.error(LeafNumber))
m1
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/LeafNumberExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Number of leaves")+
  theme_classic()+
  ylim(7,11)+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()




###### Leaf length - survey #####
fielddatS
ggplot(fielddatS, aes(x=Community,y=LeafLengthmm))+
  geom_boxplot()

m1<-fielddatS%>%
  group_by(Community)%>%
  summarise(mean=mean(LeafLengthmm,na.rm=T), se=std.error(LeafLengthmm,na.rm=T))
m1

pdf("Figs/LeafLengthSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylim(50,150)+
  ylab("Leaf length (mm)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Leaf number - Experiment #####

ggplot(fielddatE, aes(x=CommunityProj,y=LeafLengthmm))+
  geom_boxplot()

m1<-fielddatE%>%
  filter(!is.na(LeafLengthmm))%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(LeafLengthmm), se=std.error(LeafLengthmm))
m1
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/LeafLengthExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Leaf length (mm)")+
  theme_classic()+
  ylim(50,150)+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


###### Flowers length - survey #####
fielddatS
ggplot(fielddatS, aes(x=Community,y=FlowersperRosette))+
  geom_boxplot()

hist(fielddat$FlowersperRosette)
hist(fielddatS$FlowersperRosette)

#taking out samples with >6 flowers
fielddatSb<-fielddatS%>%
  filter(FlowersperRosette<6)

m1<-fielddatSb%>%
  group_by(Community)%>%
  summarise(mean=mean(FlowersperRosette,na.rm=T), se=std.error(FlowersperRosette,na.rm=T))
m1

pdf("Figs/FlowersSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylim(0,3)+
  ylab("Number of flowers")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
dev.off()



###### Flowers - Experiment #####

ggplot(fielddatE, aes(x=CommunityProj,y=FlowersperRosette))+
  geom_boxplot()

hist(fielddatE$FlowersperRosette)

#taking out samples with >6 flowers
fielddatEb<-fielddatE%>%
  filter(FlowersperRosette<6)

m1<-fielddatEb%>%
  filter(!is.na(LeafLengthmm))%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(FlowersperRosette,na.rm=T), se=std.error(FlowersperRosette,na.rm=T))
m1
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/FlowersExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Number of flowers")+
  theme_classic()+
  ylim(0,3)+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()






##### LICOR #####

licordat<-read.csv("Data/LiCordata.csv",stringsAsFactors = T)

head(licordat)

licordat2<-licordat%>%
  dplyr::select(SampleID:Community,LeafAreacm2:gtc)%>%
  filter(!is.na(LeafAreacm2))%>%
  filter(Acorrected<40)%>%
  mutate(Ecorrected=Ecorrected*1000)

licordat2$Community<-factor(licordat2$Community,levels=c("SB","WM","MM","DM","FF"))

#duplicated(licordat2$SampleID)

#TS-FF6 leaf area is small, AS-FF1 leaf area small, LS-FF1 is terrible pic, TS-FF2 has deschampsia in the chamber and is folded, EC-DM3b is folded
#filter(licordat2,Acorrected>40)

ggplot(licordat2,aes(y=Ecorrected,x=Community))+
  geom_boxplot()
ggplot(licordat2,aes(y=Acorrected,x=Community))+
  geom_boxplot()
ggplot(licordat2,aes(y=gsw,x=Community))+ #stomatal conductance
  geom_boxplot()


# Make two dataframes: for the survey data and the experiment data

licordat2S<-licordat2%>%
  filter(Proj!="E")
#filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one

licordat2E<-licordat2%>%
  filter(Proj!="S")%>%
  mutate(CommunityProj=factor(paste(Community,Proj,sep="_"),levels=c("WM_C","WM_E","MM_C","MM_E","DM_C","DM_E")))%>%
  mutate(Community=factor(Community,levels=c("WM","MM","DM")))





###### Carbon fixation - Survey ######
m1<-licordat2S%>%
  group_by(Community)%>%
  summarise(mean=mean(Acorrected), se=std.error(Acorrected))
m1

ylab(bquote('Assimilation ('*mu~'mol' ~ CO[2]~ m^-2~s^-1*')'))

pdf("Figs/PsSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  theme_classic()+
  theme(line=element_line(linewidth=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",linewidth=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  ylim(14,31)+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,linewidth=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  #ylab(bquote('Assimilation ('*mu*'mol' ~ CO[2]~ m^-2~s^-1*')'))+
  ylab(bquote("Assimilation ("*mu*"mol"~ CO[2]~m^-2~s^-1*")"))+
  geom_text(aes(y = mean+se, label = c("ab","ab","ab","a","b"),x = Community),colour="black", size=2.8,vjust = -1)
dev.off()

mod1<-lme(Acorrected~Community,random=~1|Site/Plot,weights=varIdent(form=~1|Community),data=licordat2S)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Community)))
anova(mod1, type="marginal")
summary(glht(mod1, linfct = mcp(Community = "Tukey")))

ggplot(data=mod1emm, aes(x=Community, y=emmean))+   
  geom_errorbar(aes(ymax=emmean+SE,ymin=emmean-SE),width=.2,linewidth=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab(bquote("Assimilation (\u03BCmol"~ CO[2]~m^-2~s^-1*")"))+
  theme_classic()+
  theme(line=element_line(linewidth=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",linewidth=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")




###### Carbon fixation - Experiment ######
m1<-licordat2E%>%
  #group_by(Site,Community,Proj)%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(Acorrected), se=std.error(Acorrected))
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


pdf("Figs/PsExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean,group=Community))+   
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  #ylim(14,31)+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab(bquote('Assimilation ('*mu*'mol' ~ CO[2]~ m^-2~s^-1*')'))+
  #geom_text(aes(y = mean+se, label = c("ab","ab","ab","a","b"),x = Community),colour="black", size=2.8,vjust = -1)+ 
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


mod1<-lme(Acorrected~Proj*Community,random=~1|Site/Plot,weights=varIdent(form=~1|Community),data=licordat2E)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Proj*Community)))
anova(mod1, type="marginal")
#summary(glht(mod1, linfct = mcp(Community= "Tukey")))

#the random effects are adding a lot of variance to the means
mod1<-gls(Acorrected~Community,data=licordat2E)#
summary(glht(mod1, linfct = mcp(Community= "Tukey")))



###### Transpiration - Survey ######

m1<-licordat2S%>%
  group_by(Community)%>%
  summarise(mean=mean(Ecorrected), se=std.error(Ecorrected))
m1

pdf("Figs/TranspirationSurvey.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=Community, y=mean))+   
  theme_classic()+
  theme(line=element_line(linewidth=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",linewidth=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  ylim(10,19)+
  ylab(bquote("Transpiration (mmol"~ m^-2~s^-1*")"))+
  #ylab(bquote("Transpiration (mmol"~ H[2]*"O"~m^-2~s^-1*")"))+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,linewidth=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  geom_text(aes(y = mean+se, label = c("a","abc","ab","bc","c"),x = Community),colour="black", size=2.8,vjust = -1)
dev.off()

mod1<-lme(Ecorrected~Community,random=~1|Site/Plot,weights=varIdent(form=~1|Community),data=licordat2S)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Community)))
anova(mod1, type="marginal")
summary(glht(mod1, linfct = mcp(Community = "Tukey")))




###### Transpiration - Experiment ######

m1<-licordat2E%>%
  #group_by(Site,Community,Proj)%>%
  group_by(Community,Proj)%>%
  summarise(mean=mean(Ecorrected), se=std.error(Ecorrected))
m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")

hist(licordat2E$Ecorrected)
hist(licordat2$Ecorrected)
sort(licordat2$Ecorrected)
licordat2E[licordat2E$Ecorrected>40,]
licordat2E%>%filter(Community=="MM",Proj=="E")
#there are two high points ~30 and two low points ~8 in the MM treatment

pdf("Figs/TranspirationExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Proj, y=mean,group=Community))+   
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  #ylim(14,31)+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab(bquote('Transpiration (mmol'~ m^-2~s^-1*')'))+
  #geom_text(aes(y = mean+se, label = c("ab","ab","ab","a","b"),x = Community),colour="black", size=2.8,vjust = -1)+ 
  facet_wrap(vars(Community),strip.position = "bottom")
dev.off()


mod1<-lme(Ecorrected~Proj*Community,random=~1|Site/Plot,weights=varIdent(form=~1|Community),data=licordat2E)#,correlation=corSpher(form = ~ Lat+Long)
mod1emm<-as.data.frame(summary(emmeans(mod1,~Proj*Community)))
anova(mod1, type="marginal")

mod1<-lme(Ecorrected~CommunityProj,random=~1|Site/Plot,weights=varIdent(form=~1|Community),data=licordat2E)#,correlation=corSpher(form = ~ Lat+Long)
#summary(glht(mod1, linfct = mcp(CommunityProj=c("WM_C-WM_E=0","MM_C-MM_E=0","DM_C-DM_E=0"))), test = adjusted("fdr"))
summary(glht(mod1, linfct = mcp(CommunityProj=c("WM_C-WM_E=0","MM_C-MM_E=0","DM_C-DM_E=0"))))








#### 2022 annual report ####

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

pdf("Figs/MoistureExperiment.pdf",width=2.2,height=2.2)
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

