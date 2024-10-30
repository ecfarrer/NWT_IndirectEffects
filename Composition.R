#Composition analyses


#Data files
datS
datITSS5cS
dat16SS5cS
datE
datITSS5cE
dat16SS5cE


##### Ordination Survey #####

tempphyS<-dat16SS5cS%>%
  subset_samples(SampleType=="soil")%>%
  #subset_samples(Site=="Trough")%>%
  subset_samples(CommunityType!="WM")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsS <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~CommunityType*Site))
anova(mynmdsS,by="terms",permutations = how(blocks=sample_data(tempphyS)$Site,nperm=999))
 
# mynmdsplot <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))
mynmdsplot <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))

plot_ordination(tempphyS, mynmdsplot, type="samples", color="CommunityType",axes=c(3,4))+
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

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiSurveyCommunity.pdf",width=4.3,height=3)
plot_ordination(tempphyS, mynmdsS, type="samples", color="CommunityType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityType),level=.95)
#dev.off()





##### Ordination Experiment #####
#For the experiment, Audubon is only MM, EastKnoll and Lefty are only DM, and Trough is only WM
#in elevation it goes trough, lefty, audubon, east knoll

#for 16S roots, I like the NMDS because it puts things only in 2 dimensions and especially if you plotted convex hulls you'd be able to see that there was a treatment effect in trough and audubon, but not in lefty and eastknoll. the SiteTreatment with condition(site) looks nice
#for 16S soil, all sig except east knoll. sitetreatment looks bad the sites are too different. nmds better but sites are still really different and taking up all the space. the condition on site with sitetreatment is pretty good.
#for ITS roots, the capscale with sitetreatment looks better than the nmds. all should be significant but in different directions. the SiteTreatment with condition(site) looks nice
#fot ITS soil, lefty is not sig, east knoll is p=0.045, others are more sig. doing capscale with sitetreatment or doing nmds on all you can't see an effect in trough. the SiteTreatment with condition(site) looks nice

tempphyE<-datITSS5cE%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x)) #standardizing hardly changes anything, not even noticable
sample_data(tempphyE)$Site<-factor(sample_data(tempphyE)$Site,levels=c("Trough","Audubon","Lefty","EastKnoll"))

sample_sums(tempphyE)

#mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Site*Treatment))
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~SiteTreatment+Condition(Site)))
#mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~CommunityPlotType))
anova(mynmdsE,by="terms",permutations = how(blocks=sample_data(tempphyE)$Site,nperm=999))

#mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~SiteTreatment))
#mynmdsE <- ordinate(tempphyE, "NMDS", "bray",try=200,trymax=200)

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungi.pdf",width=4.8,height=3)
plot_ordination(tempphyE, mynmdsE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)+
  facet_wrap(~Site)
#dev.off()

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungibysite.pdf",width=8,height=6)
plot_ordination(tempphyE, mynmdsE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)+
  facet_wrap(~Site)
#dev.off()



##### Experiment - separate ordinations by site, then merge #####
#Separate ordinations by site and then plot them in 4 panels. I think what is happening is that the species comps are so different among sites and the sites are all doing different things in response to the treatment, that it is impossible to show the results for all sites on the same 2 axes.

tempphyE<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")%>%
  subset_samples(Site=="Trough")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsplotE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Treatment))
anova(mynmdsplotE,by="terms",permutations = how(nperm=999))
tr<-plot_ordination(tempphyE, mynmdsplotE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+  
  #theme(legend.position = "none")+
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)
  
tempphyE<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")%>%
  subset_samples(Site=="Audubon")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsplotE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Treatment))
anova(mynmdsplotE,by="terms",permutations = how(nperm=999))
au<-plot_ordination(tempphyE, mynmdsplotE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)

tempphyE<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")%>%
  subset_samples(Site=="Lefty")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsplotE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Treatment))
anova(mynmdsplotE,by="terms",permutations = how(nperm=999))
le<-plot_ordination(tempphyE, mynmdsplotE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)

tempphyE<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")%>%
  subset_samples(Site=="EastKnoll")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsplotE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~Treatment))
anova(mynmdsplotE,by="terms",permutations = how(nperm=999))
ek<-plot_ordination(tempphyE, mynmdsplotE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)

plot_grid(tr,au,le,ek, nrow = 2,labels=c("a) Trough","b) Audubon","c) Lefty","d) East Knoll"),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain")#legend



#Community Type and treatment
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~CommunityType*Treatment))
anova(mynmdsE,by="terms",nperm=999)

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiExpCommunityType.pdf",width=8,height=3)
plot_ordination(tempphyE, mynmdsE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)+
  facet_wrap(~CommunityType)
#dev.off()



  



##### Combined survey and experiment by site #####

#If i want to do this better I need to make a single phyloseq object with all sites but with deleting the plots in each site that are necessary. it would just take a bit to do this

#In general, not surprisingly, the expermimental plots are not moving toward or overlapping the next warmer/drier community type, they are going off in their own direction

#Trough, the WM treatment plots are not moving toward the MM plots
dat16SS5cTrough
datITSS5cTrough

tempphyE<-datITSS5cTrough%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
#sample_data(tempphyE)$SampleTypeTreatment<-paste(sample_data(tempphyE)$SampleType,sample_data(tempphyE)$Treatment,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsE,by="margin",nperm=999)

plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)


#Audubon
datITSS5cAudubon

tempphyE<-dat16SS5cAudubon%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
#sample_data(tempphyE)$SampleTypeTreatment<-paste(sample_data(tempphyE)$SampleType,sample_data(tempphyE)$Treatment,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsE,by="margin",nperm=999)

plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)

#Lefty
tempphyE<-datITSS5cLefty%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
#sample_data(tempphyE)$SampleTypeTreatment<-paste(sample_data(tempphyE)$SampleType,sample_data(tempphyE)$Treatment,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsE,by="margin",nperm=999)

plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)

#East Knoll
tempphyE<-dat16SS5cEastKnoll%>%
  subset_samples(SampleType=="soil")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsE,by="margin",nperm=999)

plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)

