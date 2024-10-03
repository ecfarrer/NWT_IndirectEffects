#Composition analyses


#Data files
datITSS5
datITSS5c
datITSS5otu
datITSS5cotu

head(sample_data(datITSS5c))
sample_data(datITSS5c)$CommunityType<-factor(sample_data(datITSS5c)$CommunityType,levels=c("SB","WM","MM","DM","FF"))


sample_data(datITSS5c)$SiteCommunityType<-paste(sample_data(datITSS5c)$Site,sample_data(datITSS5c)$CommunityType,sep="")
sample_data(datITSS5c)$SiteCommunityType<-factor(sample_data(datITSS5c)$SiteCommunityType,levels=c("AudubonSB","AudubonMM","AudubonDM","AudubonFF","EastKnollSB","EastKnollMM","EastKnollDM","EastKnollFF","LeftySB","LeftyMM","LeftyDM","LeftyFF","TroughSB","TroughWM","TroughMM","TroughDM","TroughFF"))

tempphyS<-datITSS5c%>%
  subset_samples(SurveyAnalysis=="Survey")%>%
  subset_samples(SampleType=="soil")%>%
  #subset_samples(PlotType%in%c("Survey","Control"))%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

#mynmdsS <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~CommunityType*Site))

#anova(mynmdsS,by="terms",permutations = how(blocks=sample_data(tempphyS)$Site,nperm=999))

mynmdsplot <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))

plot_ordination(tempphyS, mynmdsplot, type="samples", color="CommunityType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#b30000","#ff8080"),labels = c("SB", "WM","MM","DM","FF"),name = "Community Type")+#,"#79d279"
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#b30000","#ff8080"),labels = c("SB", "WM","MM","DM","FF"),name = "Community Type")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityType),level=.95)+
  facet_wrap(~Site)

plot_ordination(tempphyS, mynmdsplot, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)


#By community type
mynmdsS <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~CommunityType+Condition(Site)))
anova(mynmdsS,by="terms",permutations = how(blocks=sample_data(tempphyS)$Site,nperm=999))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiSurveyCommunity.pdf",width=4.3,height=3)
plot_ordination(tempphyS, mynmdsS, type="samples", color="CommunityType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityType),level=.95)
dev.off()





##### Ordination Experiment #####
#For the experiment, Audubon is only MM, EastKnoll and Lefty are only DM, and Trough is only WM
#in elevation it goes trough, lefty, audubon, east knoll

head(sample_data(datITSS5c))
head(sample_data(tempphyE))

sample_data(datITSS5c)$PlotType
sample_data(datITSS5c)$Treatment

tempphyE<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
sample_data(tempphyE)$CommunityPlotType<-paste(sample_data(tempphyE)$CommunityType,sample_data(tempphyE)$PlotType,sep="")
sample_data(tempphyE)$SitePlotType<-paste(sample_data(tempphyE)$Site,sample_data(tempphyE)$PlotType,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType*Site))
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType+Condition(Site)))
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~CommunityPlotType))
anova(mynmdsE,by="terms",permutations = how(blocks=sample_data(tempphyE)$Site,nperm=999))

#mynmdsplotE <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungi.pdf",width=4.8,height=3)
plot_ordination(tempphyE, mynmdsE, type="samples", color="SitePlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=SitePlotType),level=.95)
#facet_wrap(~Site)
dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungibysite.pdf",width=8,height=6)
plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)+
  facet_wrap(~Site)
dev.off()



#Community Type and Plot Type
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~CommunityType*PlotType))
anova(mynmdsE,by="margin",nperm=999)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiExpCommunityType.pdf",width=8,height=3)
plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)+
  facet_wrap(~CommunityType)
dev.off()


#One site
tempphyE<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  #subset_samples(SampleType=="soil")%>%
  subset_samples(Site=="Audubon")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

sample_data(tempphyE)$SampleTypeTreatment<-paste(sample_data(tempphyE)$SampleType,sample_data(tempphyE)$Treatment,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Treatment+SampleType))
anova(mynmdsE,by="margin",nperm=999)

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~SampleTypeTreatment))

plot_ordination(tempphyE, mynmdsE, type="samples", color="SampleTypeTreatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=SampleTypeTreatment),level=.95)
  


