###### Plant Performance ######


#variables: soil moisture, biomass, number of leaves, leaf length, number of flowers, assimilation, transpiration, stomatal conductance (gsw)?, arbuscules, vesicles, hyphae?


options(contrasts=c("contr.treatment","contr.poly"))
options(contrasts=c("contr.helmert","contr.poly"))
options("contrasts") #see what contrasts are set

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}


###### Moisture - Survey ######

ggplot(datS, aes(x=Site,y=SoilMoisturePercent))+
  geom_boxplot()

m1<-datS%>%
  group_by(CommunityType)%>%
  summarise(mean=mean(SoilMoisturePercent), se=std.error(SoilMoisturePercent))
m1

#pdf("Figs/MoistureSurveynew.pdf",width=2.2,height=2.2)
ggplot(data=m1, aes(x=CommunityType, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  #ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")
#dev.off()

ggplot(data=datS, aes(x=CommunityType, y=SoilMoisturePercent))+   
  geom_point(data=datS,size=1.8,aes(color=CommunityType,group=CommunityType),position = position_dodge(0.8))
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Soil moisture %")+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
    #ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")

#code from mangroves for standard error bars and points
prbi<-
    ggplot(mangrove, aes(x=Substrate, y=RBI,color=Inoculum)) +
    labs(x = "Substrate",y="Root branching intensity (tips/cm)") +
    theme_classic()+
    theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
    geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
    geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
    geom_errorbar(aes(ymin = RBI-se, ymax = RBI+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
    scale_color_manual(values = c("#6c66be", "#8ca54f")) 
  #  annotate("text",x = c(1,2,3), y = 95, label = c("a","b","b"),size=3)
  
  
#Statistics

test<-datS%>%
  #filter(Site%in%c("Lefty","EastKnoll"))%>%
  filter(CommunityType!="WM")
test$CommunityType<-factor(test$CommunityType,levels=c("SB","MM","DM","FF"))

m0<-gls(SoilMoisturePercent ~ CommunityType, data = datS)
anova(m0,type="margin")

#can't have unused levels in factors
mc<-gls(SoilMoisturePercent~CommunityType,na.action=na.omit,data=datS)
summary(glht(mc, linfct = mcp(CommunityType = "Tukey")))

summary(glht(mc, linfct = mcp(SubstrateInoculum=c("GlassAutoclaved-GlassLive=0","MixAutoclaved-MixLive=0","DredgeAutoclaved-DredgeLive=0"))))
summary(glht(mc, linfct = mcp(SubstrateInoculum=c("GlassAutoclaved-DredgeAutoclaved=0","GlassAutoclaved-MixAutoclaved=0","DredgeAutoclaved-MixAutoclaved=0","GlassLive-DredgeLive=0","GlassLive-MixLive=0","DredgeLive-MixLive=0"))))



###### Moisture - Experiment #####

ggplot(datE, aes(x=MoistureTreatment,y=SoilMoisturePercent))+
  geom_boxplot()

m1<-datE%>%
  group_by(MoistureType,Treatment)%>%
  summarise(mean=mean(SoilMoisturePercent), se=std.error(SoilMoisturePercent))
m1
#m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


#pdf("Figs/MoistureExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Treatment, y=mean))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+
  ylab("Soil moisture %")+
  #ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(MoistureType),strip.position = "bottom")
#dev.off()


#Statistics

m0<-gls(SoilMoisturePercent ~ Treatment*MoistureType,  data = datE)
anova(m0,type="margin")

m0<-lme(SoilMoisturePercent ~ Treatment, random=~1|Site/Treatment/Plot, data = datE)
anova(m0,type="margin")





###### Biomass ######

ggplot(datE, aes(x=MoistureTreatment,y=Biomassg))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(Biomassg),Biomassg=mean(Biomassg))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=Biomassg, color=Treatment))+   
  ylab("Biomass (g)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=Biomassg+se,ymin=Biomassg-se),width=.3)+
  ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
#  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(Biomassg ~ Treatment*CommunityType, random=~1|SitePlot/Chamber, data = datE)
anova(m0,type="margin")
mc<-lme(Biomassg ~ SiteTreatment, random=~1|SitePlot/Chamber, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))


###### Leaf number ######

ggplot(datE, aes(x=Treatment,y=LeafNumber))+
  geom_boxplot()
hist(datE$LeafNumber)

datE$LeafNumber[datE$LeafNumber>20]<-NA ##i replaced the 21 with NA

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(LeafNumber),LeafNumber=mean(LeafNumber,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=LeafNumber, color=Treatment))+   
  ylab("Leaf number")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=LeafNumber+se,ymin=LeafNumber-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(LeafNumber ~ Treatment*Site, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(LeafNumber ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))


###### Leaf length ######

ggplot(datE, aes(x=Treatment,y=LeafLengthmm))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(LeafLengthmm),LeafLengthmm=mean(LeafLengthmm))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=LeafLengthmm, color=Treatment))+   
  ylab("Leaf number")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=LeafLengthmm+se,ymin=LeafLengthmm-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(LeafLengthmm ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(LeafLengthmm ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))



###### FlowersperRosette ######

ggplot(datE, aes(x=Treatment,y=FlowersperRosette))+
  geom_boxplot()
hist(datE$FlowersperRosette)

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(FlowersperRosette),FlowersperRosette=mean(FlowersperRosette,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=FlowersperRosette, color=Treatment))+   
  ylab("Flowers per rossette")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=FlowersperRosette+se,ymin=FlowersperRosette-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(FlowersperRosette ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(FlowersperRosette ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))



###### Assimilation ######

ggplot(datE, aes(x=Treatment,y=Acorrectedredo))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(Acorrectedredo),Acorrectedredo=mean(Acorrectedredo,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=Acorrectedredo, color=Treatment))+   
  ylab(bquote('Assimilation ('*mu*'mol' ~ CO[2]~ m^-2~s^-1*')'))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=Acorrectedredo+se,ymin=Acorrectedredo-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(Acorrectedredo ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(Acorrectedredo ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))


###### Transpiration ######

ggplot(datE, aes(x=Treatment,y=Ecorrectedredo))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(Ecorrectedredo),Ecorrectedredo=mean(Ecorrectedredo,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=Ecorrectedredo, color=Treatment))+   
  ylab(bquote('Transpiration (mmol'~ m^-2~s^-1*')'))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=Ecorrectedredo+se,ymin=Ecorrectedredo-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(Ecorrectedredo ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(Ecorrectedredo ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))


###### Stomatal conductance ######

ggplot(datE, aes(x=Treatment,y=gsw))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(gsw),gsw=mean(gsw,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=gsw, color=Treatment))+   
  ylab(bquote('Stomatal conductance (mol'~ m^-2~s^-1*')'))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=gsw+se,ymin=gsw-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(gsw ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(gsw ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))



###### gtc (total conductance to CO2), gtw (total conductance to water vapor), Ci (intracellular CO2) ######

#The conductances are all pretty much the same exact pattern as gsw. Ci is interesting,there is a site effect with trough high and audubon low. a nonsig positive effect of treatment at Audubon
m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(Ci),Ci=mean(Ci,na.rm=T))
m1
#m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=Ci, color=Treatment))+   
  ylab(bquote('Stomatal conductance (mol'~ m^-2~s^-1*')'))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.25),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=Ci+se,ymin=Ci-se),width=.3)+
  #ylim(0,2.8)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  #  geom_text(m1,mapping=aes(x=Treatment,y=2.7,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(Ci ~ Site*Treatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
anova(m0,type="margin")
mc<-lme(Ci ~ SiteTreatment, random=~1|SitePlot/Chamber, na.action=na.omit, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))



###### Arbuscules ######

ggplot(datE, aes(x=MoistureTreatment,y=ArbusculesP))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise( se=std.error(ArbusculesP),ArbusculesP=mean(ArbusculesP))
m1
m1$mcs<-c("a","b","a","a","a","a","a","a")

#with facet wrap
ggplot(datE, aes(x=Treatment, y=ArbusculesP, color=Treatment))+   
  ylab("Arbuscules %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=ArbusculesP+se,ymin=ArbusculesP-se),width=.3,position=position_dodge(0.8))+
  ylim(0,110)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  geom_text(m1,mapping=aes(x=Treatment,y=105,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-lme(ArbusculesP ~ Treatment*Site, random=~1|SitePlot/Chamber, data = datE)
anova(m0,type="margin")
mc<-lme(ArbusculesP ~ SiteTreatment, random=~1|SitePlot/Chamber, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))



###### Vesicles ######

ggplot(datE, aes(x=MoistureTreatment,y=VesiclesP))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise( se=std.error(VesiclesP),VesiclesP=mean(VesiclesP))
m1
m1$mcs<-c("a","b","a","a","a","a","a","a")
#m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")
#I could convert Trough=1, Audubon=2,Lefty=3, EastKnoll=4 and then use scale_x_continuous to mess with labels. I just have to be careful with relabeling the sites correctly

#no facet wrap, working but can't get treatment on the x axis
# ggplot(datE, aes(x=Site, y=VesiclesP, color=Treatment))+   
#   ylab("Vesicles %")+
#   theme_classic()+
#   theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "right",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
#   geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment),position=position_dodge(0.8))+
#   geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,72))+
#   geom_errorbar(data=m1,aes(ymax=VesiclesP+se,ymin=VesiclesP-se),width=.3,position=position_dodge(0.8))+
#   annotate("text",x = c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2), y = 95, label = c("a","b","a","a","a","a","a","a"),size=3)+  
#   scale_color_manual(values = c("#6c66be", "#8ca54f")) 

#with facet wrap
ggplot(datE, aes(x=Treatment, y=VesiclesP, color=Treatment))+   
  ylab("Vesicles %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=VesiclesP+se,ymin=VesiclesP-se),width=.3,position=position_dodge(0.8))+
  ylim(0,110)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  geom_text(m1,mapping=aes(x=Treatment,y=105,label=mcs),color="black",size=3)+
  scale_color_manual(values = c("#50b47b", "#ba6437")) 


#Statistics
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-gls(VesiclesP ~ Treatment*MoistureType,  data = datE)
anova(m0,type="margin")

m0<-lme(VesiclesP ~ Treatment*MoistureType, random=~1|Site/Treatment/Plot, data = datE)
anova(m0,type="margin")

m0<-lme(VesiclesP ~ Treatment*MoistureType, random=~1|Site/Treatment, data = datE)
anova(m0,type="margin")

m0<-lme(VesiclesP ~ Treatment*Site, random=~1|SitePlot/Chamber, data = datE)
anova(m0,type="margin")
mc<-lme(VesiclesP ~ SiteTreatment, random=~1|SitePlot/Chamber, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))


###### Vesicles/Arbuscules from undergrads ######
#upshot: the patterns in arbuscules and vesicles are kind of similar to Monica's. the undergrads see a slight increase in Trough vesicles with warming (but not nearly as much as Monica/anna's data). they see a slighly larger increase in Audubon vesicles with warming compared to Monica's data. This suggests I could average the vesicle data for Trough and Audubon to make it more robust to multiple scorers. However the arbuscule data is much lower in the undergrad dataset possibly the stain was leeching out or they had a hard time seeing the arbuscules (overall average of 15% vs 50%). the pattern of more arbuscules in the warmed plots is shown in both the undergrad and  Monica's dataset for Tough and Audubon. but again the difference is greater in Monica's trough vs undergrad trough, suggesting that the Anna/Monica scorer thing is a problem. the difference in audubon is very similar between Monica vs undergrads.

uamf2
ggplot(uamf2, aes(x=SiteTreatment,y=ArbusculesP))+
  geom_boxplot()

m1<-uamf2%>%
  filter(ExperimentAnalysis=="Experiment")%>%
  filter(Site!="EastKnoll")%>%
  group_by(Site,Treatment)%>%
  summarise(se=std.error(ArbusculesP),ArbusculesP=mean(ArbusculesP))
m1
m1$mcs<-c("a","a","a","a","a","a")

#with facet wrap
ggplot(uamf2, aes(x=Treatment, y=ArbusculesP, color=Treatment))+   
  ylab("Vesicles %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment))+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),alpha=rep(.4,68))+
  geom_errorbar(data=m1,aes(ymax=ArbusculesP+se,ymin=ArbusculesP-se),width=.3,position=position_dodge(0.8))+
  #ylim(0,110)+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)#+
#  geom_text(m1,mapping=aes(x=Treatment,y=105,label=mcs),color="black",size=3)+
#  scale_color_manual(values = c("#50b47b", "#ba6437")) 



