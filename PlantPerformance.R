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




###### Vesicles ######

ggplot(datE, aes(x=MoistureTreatment,y=VesiclesP))+
  geom_boxplot()

m1<-datE%>%
  group_by(Site,Treatment)%>%
  summarise( se=std.error(VesiclesP),VesiclesP=mean(VesiclesP))
m1
#m1$Proj<-recode_factor(m1$Proj,"C"="Control","E"="Treatment")


#pdf("Figs/MoistureExperiment.pdf",width=3.2,height=2.2)
ggplot(datE, aes(x=Site, y=VesiclesP, color=Treatment))+   
  ylab("Vesicles %")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  geom_point(data=m1,size=1.8,aes(color=Treatment,group=Treatment),position=position_dodge(0.8),show.legend = FALSE)+
  geom_point(size=.7,aes(fill=Treatment,group=Treatment),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,72))+
  geom_errorbar(data=m1,aes(ymax=VesiclesP+se,ymin=VesiclesP-se),width=.3,position=position_dodge(0.8))+
  #ylim(5,60)+
  #facet_wrap(vars(Site),strip.position = "bottom",nrow=1)+
  annotate("text",x = c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2), y = 70, label = c("a","b","a","a","a","a","a","a"),size=3)+  
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
#dev.off()

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
#Chamber is unique to the 1 x 1 m plot, "Plot" is unique to the elevation of the site across sand treatments (1,2,3 1,2,3)

m0<-gls(Vesicles ~ Treatment*MoistureType,  data = datE)
anova(m0,type="margin")

m0<-lme(Vesicles ~ Treatment*MoistureType, random=~1|Site/Treatment/Plot, data = datE)
anova(m0,type="margin")

m0<-lme(Vesicles ~ Treatment*MoistureType, random=~1|Site/Treatment, data = datE)
anova(m0,type="margin")

m0<-lme(Vesicles ~ Treatment*Site, random=~1|SitePlot/Chamber, data = datE)
anova(m0,type="margin")
mc<-lme(Vesicles ~ SiteTreatment, random=~1|SitePlot/Chamber, data = datE)
summary(glht(mc, linfct = mcp(SiteTreatment=c("Trough_Experimental-Trough_Control=0","Audubon_Experimental-Audubon_Control=0","Lefty_Experimental-Lefty_Control=0","EastKnoll_Experimental-EastKnoll_Control=0"))))

