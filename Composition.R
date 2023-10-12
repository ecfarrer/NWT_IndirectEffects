#Composition analyses

options(contrasts=c("contr.treatment","contr.poly"));options("contrasts")
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

#Data files
datITSS5
datITSS5c
datITSS5otu
datITSS5cotu

head(sample_data(datITSS5c))
sample_data(datITSS5c)$CommunityType<-factor(sample_data(datITSS5c)$CommunityType,levels=c("SB","WM","MM","DM","FF"))

#### Ordination ####


##### Ordination survey #####

head(sample_data(datITSS5c))

sample_data(datITSS5c)$SiteCommunityType<-paste(sample_data(datITSS5c)$Site,sample_data(datITSS5c)$CommunityType,sep="")
sample_data(datITSS5c)$SiteCommunityType<-factor(sample_data(datITSS5c)$SiteCommunityType,levels=c("AudubonSB","AudubonMM","AudubonDM","AudubonFF","EastKnollSB","EastKnollMM","EastKnollDM","EastKnollFF","LeftySB","LeftyMM","LeftyDM","LeftyFF","TroughSB","TroughWM","TroughMM","TroughDM","TroughFF"))

tempphyS<-datITSS5c%>%
  #subset_samples(SurveyAnalysis=="Survey")%>%
  subset_samples(PlotType%in%c("Survey","Control"))%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsS <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~CommunityType*Site))
anova(mynmdsS,by="terms",permutations = how(blocks=sample_data(tempphyS)$Site,nperm=999))

mynmdsplot <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/PCoAfungi.pdf",width=4.8,height=3)
plot_ordination(tempphyS, mynmdsplot, type="samples", color="CommunityType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#b30000","#ff8080"),labels = c("SB", "WM","MM","DM","FF"),name = "Community Type")+#,"#79d279"
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#b30000","#ff8080"),labels = c("SB", "WM","MM","DM","FF"),name = "Community Type")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityType),level=.95)+
  facet_wrap(~Site)
#dev.off()

plot_ordination(tempphyS, mynmdsplot, type="samples", color="CommunityType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityType),level=.95)


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

head(sample_data(datITSS5c))
head(sample_data(tempphyE))

tempphyE<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
sample_data(tempphyE)$CommunityPlotType<-paste(sample_data(tempphyE)$CommunityType,sample_data(tempphyE)$PlotType,sep="")

mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType*Site))
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~PlotType+Condition(Site)))
mynmdsE <- ordinate(tempphyE, "CAP", "bray",formula=as.formula(~CommunityPlotType))
anova(mynmdsE,by="terms",permutations = how(blocks=sample_data(tempphyE)$Site,nperm=999))

#mynmdsplotE <- ordinate(tempphyS, "CAP", "bray",formula=as.formula(~SiteCommunityType))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungi.pdf",width=4.8,height=3)
plot_ordination(tempphyE, mynmdsE, type="samples", color="CommunityPlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=CommunityPlotType),level=.95)
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




###### NMDS outputs for path analysis ######
tempphyEWM<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  subset_samples(CommunityType=="WM")%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T) #note that filtering >2 can take out as many has 3800 reads (from 7198)
#sample_data(tempphyEWM)$CommunityPlotType<-paste(sample_data(tempphyEWM)$CommunityType,sample_data(tempphyEWM)$PlotType,sep="")

rowSums(otu_table(tempphyEWM))

mynmdsE <- ordinate(tempphyEWM, "NMDS", "bray")
#mynmdsE <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~PlotType))
#anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)

scores(mynmdsE)$sites
sample_data(tempphyEWM)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWM)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEWM, "NMDS", "bray")
sample_data(tempphyEWM)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWM)$NMDS2<-scores(mynmdsE)$sites[,2]


boxplot(sample_data(tempphyEWM)$MDS2~sample_data(tempphyEWM)$PlotType)

sample_data(tempphyEWM)


tempphyEDM<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  subset_samples(CommunityType=="DM")%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T) 
#sample_data(tempphyEWM)$CommunityPlotType<-paste(sample_data(tempphyEWM)$CommunityType,sample_data(tempphyEWM)$PlotType,sep="")

rowSums(otu_table(tempphyEDM))

mynmdsE <- ordinate(tempphyEDM, "NMDS", "bray")
#mynmdsE <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~PlotType))
#anova(mynmdsE)
mynmdsE <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyE, mynmdsE, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)

scores(mynmdsE)$sites
sample_data(tempphyEDM)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEDM)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEDM, "NMDS", "bray")
sample_data(tempphyEDM)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEDM)$NMDS2<-scores(mynmdsE)$sites[,2]


boxplot(sample_data(tempphyEDM)$MDS2~sample_data(tempphyEDM)$PlotType)

sample_data(tempphyEDM)









###### Ordination of stress and abiotic affecting microbes ######

#goal to test if moisture, Ps, Transpiration, Vesicles affects microbes

tempphyE<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")
temp<-sample_data(tempphyE)
temp$SampleID<-gsub("_", "-", temp$PlotID)
temp$SampleID<-sub("(-.*?)-", "\\1", temp$SampleID)
head(temp)

temp2<-data.frame(temp)%>%
  left_join(dat,by="SampleID")
head(temp2)
dim(temp2)

#they are in the same order
temp$SampleNumber
temp2$SampleNumber

#this is wierd, the first time i did it was 18:57
#sample_data(tempphyE)[,18:57]<-temp2[,18:57]
sample_data(tempphyE)[,17:55]<-temp2[,18:56]


# WM
tempphyEWM<-tempphyE%>%
  subset_samples(Community=="WM")%>%
  subset_samples(!is.na(Acorrected))%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)
#sample_data(tempphyEMM)$CommunityPlotType<-paste(sample_data(tempphyE)$CommunityType,sample_data(tempphyE)$PlotType,sep="")

sample_data(tempphyEWM)$Acorrected

mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsEWM,by="margin",nperm=999)
mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEWM,by="margin",nperm=999)

#Backwards selection
mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEWM,by="margin",nperm=999)
mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEWM,by="margin",nperm=999)
mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~Acorrected+VesiclesP))
anova(mynmdsEWM,by="margin",nperm=999)
anova(mynmdsEWM,nperm=999)

#WM communities are explained by Ps and vesicles

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiExpCommunityType.pdf",width=8,height=3)
plot_ordination(tempphyEWM, mynmdsEWM, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)
dev.off()

plot(mynmdsEWM)
mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Condition(Acorrected+Ecorrected+VesiclesP)))
anova(mynmdsEWM,by="margin",nperm=999)

mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~Ecorrected+Condition(Acorrected+SoilMoisturePercent+VesiclesP)))
anova(mynmdsEWM,by="margin",nperm=999)

mynmdsEWM <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~VesiclesP+Acorrected+Condition(Ecorrected+SoilMoisturePercent)))
anova(mynmdsEWM,nperm=999)



# MM
tempphyEMM<-tempphyE%>%
  subset_samples(Community=="MM")%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)
#sample_data(tempphyEMM)$CommunityPlotType<-paste(sample_data(tempphyE)$CommunityType,sample_data(tempphyE)$PlotType,sep="")

sample_data(tempphyEMM)

mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsEMM,by="margin",nperm=999)
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEMM,by="margin",nperm=999)

#Backwards selection
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEMM,by="margin",nperm=999)
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEMM,by="margin",nperm=999)
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~Acorrected+Ecorrected))
anova(mynmdsEMM,by="margin",nperm=999)
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~Acorrected))
anova(mynmdsEMM,by="margin",nperm=999)

#nothing is significant

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiExpCommunityType.pdf",width=8,height=3)
plot_ordination(tempphyEMM, mynmdsEMM, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)
dev.off()

plot(mynmdsEMM)
mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Condition(Acorrected+Ecorrected+VesiclesP)))
anova(mynmdsEMM,by="margin",nperm=999)

mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~Ecorrected+Condition(Acorrected+SoilMoisturePercent+VesiclesP)))
anova(mynmdsEMM,by="margin",nperm=999)

mynmdsEMM <- ordinate(tempphyEMM, "CAP", "bray",formula=as.formula(~VesiclesP+Acorrected+Condition(Ecorrected+SoilMoisturePercent)))
anova(mynmdsEMM,nperm=999)




# DM
tempphyEDM<-tempphyE%>%
  subset_samples(Community=="DM")%>%
  subset_samples(!is.na(Acorrected))%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)
#sample_data(tempphyEMM)$CommunityPlotType<-paste(sample_data(tempphyE)$CommunityType,sample_data(tempphyE)$PlotType,sep="")

sample_data(tempphyEDM)$VesiclesP

mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~PlotType))
anova(mynmdsEDM,by="margin",nperm=999)
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEDM,by="margin",nperm=999)
anova(mynmdsEDM,nperm=999)

#Backwards selection
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Acorrected+Ecorrected+VesiclesP))
anova(mynmdsEDM,by="margin",nperm=999)
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Ecorrected+VesiclesP))
anova(mynmdsEDM,by="margin",nperm=999)
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Ecorrected))
anova(mynmdsEDM,by="margin",nperm=999)
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~Ecorrected))
anova(mynmdsEDM,by="margin",nperm=999)

#the only variable that affects DM community comp is transpiration rate


pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Figs/dbRDArootsfungiExpCommunityType.pdf",width=8,height=3)
plot_ordination(tempphyEDM, mynmdsEDM, type="samples", color="PlotType",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=PlotType),level=.95)
dev.off()

plot(mynmdsEDM)
mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~SoilMoisturePercent+Condition(Acorrected+Ecorrected+VesiclesP)))
anova((mynmdsEDM),by="margin",nperm=999)

mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~Ecorrected+Condition(Acorrected+SoilMoisturePercent+VesiclesP)))
anova(mynmdsEDM,by="margin",nperm=999)

mynmdsEDM <- ordinate(tempphyEDM, "CAP", "bray",formula=as.formula(~VesiclesP+Acorrected+Condition(Ecorrected+SoilMoisturePercent)))
anova(mynmdsEDM,nperm=999)























##### PLANTS ordination #####

#Setting up a phyloseq data set for plants

#without phrag
otuplants<-dat17[,6:51]/rowSums(dat17[,6:51])
colnames(otuplants) <- gsub(" ",".",colnames(otuplants))
rownames(otuplants)<-paste("s",dat17$Plot,sep="")
ind<-which(is.na(rowSums(otuplants)))
otuplants[ind,]<-0
sampleplants<-dat17[,c(1:4,52:81)]
rownames(sampleplants)<-paste("s",sampleplants$Plot,sep="")
datPlants <- merge_phyloseq(otu_table(otuplants, taxa_are_rows=FALSE), 
                            sample_data(sampleplants))


#with phrag, use this
otuplants<-dat17[,5:51]/rowSums(dat17[,5:51])
colnames(otuplants) <- gsub(" ",".",colnames(otuplants))
rownames(otuplants)<-paste("s",dat17$Plot,sep="")
ind<-which(is.na(rowSums(otuplants)))
otuplants[ind,]<-0
sampleplants<-dat17[,c(1:4,52:116)]
rownames(sampleplants)<-paste("s",sampleplants$Plot,sep="")

datPlants <- merge_phyloseq(otu_table(otuplants, taxa_are_rows=FALSE), 
                            sample_data(sampleplants))

#with tax table
otuplants<-dat17[,5:51]/rowSums(dat17[,5:51])
colnames(otuplants) <- gsub(" ",".",colnames(otuplants))
rownames(otuplants)<-paste("s",dat17$Plot,sep="")
ind<-which(is.na(rowSums(otuplants)))
otuplants[ind,]<-0
otuplants2<-t(otuplants)
sampleplants<-dat17[,c(1:4,52:81)]
rownames(sampleplants)<-paste("s",sampleplants$Plot,sep="")
taxonomyplants<-as.matrix(data.frame(Kingdom=row.names(otuplants2),Phylum=row.names(otuplants2),Class=row.names(otuplants2),Order=row.names(otuplants2),Class=row.names(otuplants2),Family=row.names(otuplants2),Genus=row.names(otuplants2),Species=row.names(otuplants2)))
rownames(taxonomyplants)<-row.names(otuplants2)
datPlants <- merge_phyloseq(otu_table(otuplants2,taxa_are_rows = T), tax_table(taxonomyplants), sample_data(sampleplants))


sample_data(datPlants)$Transect<-factor(sample_data(datPlants)$Transect,levels = c("Native","Transition","Phragmites"))
sites<-c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1")
ordlistp<-list(NA)
percentexplained<-data.frame(percentplant=rep(NA,8),fplant=rep(NA,8),sigplant=rep(NA,8),percentfungi=rep(NA,8),ffungi=rep(NA,8),sigfungi=rep(NA,8),percentbac=rep(NA,8),fbac=rep(NA,8),sigbac=rep(NA,8))

for (i in 1:8){
  tempphy<-datPlants%>%
    subset_samples(Site==sites[i])%>%
    filter_taxa(function(x) sum(x>0) >0, prune=T)
  tempphy<-prune_samples(sample_sums(tempphy)>0,tempphy)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistp[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,1]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,2]<-tempanova$F[1]
  percentexplained[i,3]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistp[[1]],ordlistp[[2]],ordlistp[[3]],ordlistp[[4]],ordlistp[[5]],ordlistp[[6]],ordlistp[[7]],ordlistp[[8]],legend, nrow = 1)

#Thoughts: it is the ellipse code that gives the convergence failure warning. the ellipses look pretty horrible anyway, should make them some other way (they seem not to be an ellipse around your points but rather an ellipse around the centroid with error). With almost all site you have problems with dispersion being different among plots. The BB is odd b/c one native plot is 100% junroe.


#Final analyses and figure for manuscript
tempphyP<-datPlants%>%
  subset_samples(Transect!="Transition")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
#tempphyP<-prune_samples(sample_sums(tempphyP)>0,tempphyP)
#mynmdsP <- ordinate(tempphyP, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
#mynmdsP <- ordinate(tempphyP, "CAP","bray",formula=as.formula(~MarshClass_Transect))

mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV*Transect))
#mynmdsP <- ordinate(tempphyP, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
#mynmds <- ordinate(tempphyP, "NMDS", "bray")
#mynmds <- ordinate(tempphyP, "NMDS", distance(tempphyP, method = "jaccard", binary = TRUE))
anova(mynmdsP,by="terms",permutations = how(blocks=sample_data(tempphyP)$Site,nperm=9999))
mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV+Transect))
anova(mynmdsP,by="margin",permutations = how(blocks=sample_data(tempphyP)$Site,nperm=9999)) #permuting within site doesn't change the p value or f statistic

#Variance explained
mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsP
mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsP
mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV*Transect+Condition(Transect+MarshClassV)));mynmdsP

mynmdsP <- ordinate(tempphyP, "CAP",distance(tempphyP, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassVTransect))
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/PCoAplant.pdf",width=4.8,height=3)
plot_ordination(tempphyP, mynmdsP, type="samples", color="MarshClassVTransect",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #geom_mark_ellipse(expand = 0,aes(fill=MarshClass_Transect))
  scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=MarshClass_Transect),level=.95)
dev.off()

#fresh 0066ff
#brackish 669900
#saline cc0000


##### Only on haplotype I ####
tempphyPI<-tempphyP%>%
  subset_samples(Site!="LUMCON 1",Site!="LUMCON 2")
mynmdsPI <- ordinate(tempphyPI, "CAP",distance(tempphyPI, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV*Transect))
anova(mynmdsPI,by="terms",permutations = how(blocks=sample_data(tempphyPI)$Site,nperm=9999))
mynmdsPI <- ordinate(tempphyPI, "CAP",distance(tempphyPI, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV+Transect))
anova(mynmdsPI,by="margin",permutations = how(blocks=sample_data(tempphyPI)$Site,nperm=9999))
mynmdsPI <- ordinate(tempphyPI, "CAP",distance(tempphyPI, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsPI
mynmdsPI <- ordinate(tempphyPI, "CAP",distance(tempphyPI, method = "jaccard", binary = TRUE),formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsPI
mynmdsPI <- ordinate(tempphyPI, "CAP",distance(tempphyPI, method = "jaccard", binary = TRUE),formula=as.formula(~MarshClassV*Transect+Condition(Transect+MarshClassV)));mynmdsPI

#Testing out varpart
tempphyPIdf<-data.frame(sample_data(tempphyPI))
#varpart still does not allow you to test the marginal effect of a main effect over the interaction effect
varpart(distance(tempphyPI, method = "jaccard", binary = TRUE), ~MarshClassV,~Transect,data=tempphyPIdf)
varpart(t(otu_table(tempphyPI)), ~MarshClassV,~Transect,~MarshClassV*Transect,data=tempphyPIdf)


#upshot, need to use with phrag, either bray and nmds or jaccard and pcoa


##### Plants - separate ordinations with jaccard including phrag by each marsh type ####

datPlantsF<-datPlants%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Fresh")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)

mynmdsFresh <- ordinate(datPlantsF, "CAP",distance(datPlantsF, method = "jaccard", binary = TRUE),formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsFresh,permutations = how(blocks=sample_data(datPlantsF)$Site,nperm=9999),by="margin")
mynmdsFresh <- ordinate(datPlantsF, "CAP",distance(datPlantsF, method = "jaccard", binary = TRUE),formula=as.formula(~1+Condition(Site)))#
Plant1<-plot_ordination(datPlantsF, mynmdsFresh, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  theme(legend.position = "none")+
  geom_point(size = 2)+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

sppscores(mynmdsFresh) <- data.frame(t(otu_table(datPlantsF)))
plot(mynmdsFresh)
ind<-order(scores(mynmdsFresh, display = "species")[,1])
scores(mynmdsFresh, display = "species")[ind,]


datPlantsB<-datPlants%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Brackish")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)

mynmdsBrackish <- ordinate(datPlantsB, "CAP",distance(datPlantsB, method = "jaccard", binary = TRUE),formula=as.formula(~Transect+Condition(Site)))
anova(mynmdsBrackish,permutations = how(blocks=sample_data(datPlantsB)$Site,nperm=9999),by="margin")
mynmdsBrackish <- ordinate(datPlantsB, "CAP",distance(datPlantsB, method = "jaccard", binary = TRUE),formula=as.formula(~1+Condition(Site)))
Plant2<-plot_ordination(datPlantsB, mynmdsBrackish, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  theme(legend.position = "none")+
  geom_point(size = 2)+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

sppscores(mynmdsBrackish) <- data.frame(t(otu_table(datPlantsB)))
plot(mynmdsBrackish)
ind<-order(scores(mynmdsBrackish, display = "species")[,1])
scores(mynmdsBrackish, display = "species")[ind,]


datPlantsS<-datPlants%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Saline")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)

mynmdsSaline <- ordinate(datPlantsS, "CAP",distance(datPlantsS, method = "jaccard", binary = TRUE),formula=as.formula(~Transect+Condition(Site)))
anova(mynmdsSaline,permutations = how(blocks=sample_data(datPlantsS)$Site,nperm=9999),by="margin")
mynmdsSaline <- ordinate(datPlantsS, "CAP",distance(datPlantsS, method = "jaccard", binary = TRUE),formula=as.formula(~1+Condition(Site)))
Plant3<-plot_ordination(datPlantsS, mynmdsSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  theme(legend.position = "none")+
  geom_point(size = 2)+  
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

legendPlant<-plot_ordination(datPlantsS, mynmdsSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+  
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/SeparateOrdinationsPlant.pdf",width=12.4,height=3)#4.8
legend <- cowplot::get_legend(legendPlant)
plot_grid(Plant1,Plant2,Plant3,legend,nrow = 1)
dev.off()


##### FUNGI #####
ordlistf<-list(NA)
for (i in 1:8){
  tempphy<-datITSS5c%>%
    subset_samples(Site==sites[i])%>%
    filter_taxa(function(x) sum(x>0) >2, prune=T)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistf[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,4]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,5]<-tempanova$F[1]
  percentexplained[i,6]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistf[[1]],ordlistf[[2]],ordlistf[[3]],ordlistf[[4]],ordlistf[[5]],ordlistf[[6]],ordlistf[[7]],ordlistf[[8]],legend, nrow = 1)

#Trying to do a single ordination. this doesn't work b/c you can't condition on site and then have marshclassV as an explanatory variable. can do it if you ignore site effects.However, in this analysis, the sig marshclass by transition interaction means that marsh composition is changing in 

sample_data(datITSS5c)$"MarshClass_Transect"<-factor(sample_data(datITSS5c)$"MarshClass_Transect",levels=c("Fresh Native","Fresh Transition","Fresh Monoculture","Brackish Native","Brackish Transition","Brackish Monoculture","Saline Native","Saline Transition","Saline Monoculture"))

sample_data(datITSS5c)$MarshClassVTransect<-factor(sample_data(datITSS5c)$MarshClassVTransect,levels=c("Fresh Native","Fresh Transition","Fresh Phragmites","Brackish Native","Brackish Transition","Brackish Phragmites","Saline Native","Saline Transition","Saline Phragmites"))

tempphyF<-datITSS5c%>%
  subset_samples(Transect!="Transition")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
anova(mynmdsF,by="terms",permutations = how(blocks=sample_data(tempphyF)$Site,nperm=9999))
mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~MarshClassV+Transect))
anova(mynmdsF,by="margin",permutations = how(blocks=sample_data(tempphyF)$Site,nperm=9999))

#variance explained
mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsF
mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsF
mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~MarshClassV*Transect+Condition(MarshClassV+Transect)));mynmdsF


mynmdsF <- ordinate(tempphyF, "CAP", "bray",formula=as.formula(~MarshClassVTransect))
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/PCoAfungi.pdf",width=4.8,height=3)
plot_ordination(tempphyF, mynmdsF, type="samples", color="MarshClassVTransect",axes=c(1,3))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Monoculture","Brackish Native","Brackish Monoculture","Saline Native","Saline Monoculture"),name = "Marsh Class/Transect")+
  scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Monoculture","Brackish Native","Brackish Monoculture","Saline Native","Saline Monoculture"),name = "Marsh Class/Transect")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=MarshClass_Transect),level=.95)
dev.off()


##### Only on haplotype I #####
tempphyFI<-tempphyF%>%
  subset_samples(Site!="LUMCON 1",Site!="LUMCON 2")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsFI <- ordinate(tempphyFI, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
anova(mynmdsFI,by="terms",permutations = how(blocks=sample_data(tempphyFI)$Site,nperm=9999))
mynmdsFI <- ordinate(tempphyFI, "CAP", "bray",formula=as.formula(~MarshClassV+Transect))
anova(mynmdsFI,by="margin",permutations = how(blocks=sample_data(tempphyFI)$Site,nperm=9999))
mynmdsFI <- ordinate(tempphyFI, "CAP","bray",formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsFI
mynmdsFI <- ordinate(tempphyFI, "CAP","bray",formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsFI
mynmdsFI <- ordinate(tempphyFI, "CAP","bray",formula=as.formula(~MarshClassV*Transect+Condition(Transect+MarshClassV)));mynmdsFI




###### Fungi - separate ordinations within marsh type #####
tempphyFF<-datITSS5c%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Fresh")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsFFresh <- ordinate(tempphyFF, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsFFresh,permutations = how(blocks=sample_data(tempphyFF)$Site,nperm=9999))
mynmdsFFresh <- ordinate(tempphyFF, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Fungi1<-plot_ordination(tempphyFF, mynmdsFFresh, type=c("sites"),color="Transect",axes=c(1,2))+#, color="Class"
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]



tempphyFB<-datITSS5c%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Brackish")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsFBrackish <- ordinate(tempphyFB, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsFBrackish,permutations = how(blocks=sample_data(tempphyFB)$Site,nperm=9999))
mynmdsFBrackish <- ordinate(tempphyFB, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Fungi2<-plot_ordination(tempphyFB, mynmdsFBrackish, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]



tempphyFS<-datITSS5c%>%
  #subset_samples(Plot!=168)%>%  #doesn't change much
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Saline")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsFSaline <- ordinate(tempphyFS, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsFSaline,permutations = how(blocks=sample_data(tempphyFS)$Site,nperm=9999))
mynmdsFSaline <- ordinate(tempphyFS, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Fungi3<-plot_ordination(tempphyFS, mynmdsFSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]

legendPlant<-plot_ordination(datPlantsS, mynmdsSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+  
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/SeparateOrdinationsFungi.pdf",width=12.4,height=3)#4.8
legend <- cowplot::get_legend(legendPlant)
plot_grid(Fungi1,Fungi2,Fungi3,legend,nrow = 1)
dev.off()





###### BACTERIA ordination #####
ordlistb<-list(NA)
for (i in 1:8){
  tempphy<-datBac4c%>%
    subset_samples(Plot!=168)%>% #take out crazy 2 million read plot
    subset_samples(Site==sites[i])%>%
    filter_taxa(function(x) sum(x>0) >2, prune=T)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistb[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,7]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,8]<-tempanova$F[1]
  percentexplained[i,9]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistb[[1]],ordlistb[[2]],ordlistb[[3]],ordlistb[[4]],ordlistb[[5]],ordlistb[[6]],ordlistb[[7]],ordlistb[[8]],legend, nrow = 1)

sample_data(datBac4c)$"MarshClass_Transect"<-factor(sample_data(datBac4c)$"MarshClass_Transect",levels=c("Fresh Native","Fresh Transition","Fresh Monoculture","Brackish Native","Brackish Transition","Brackish Monoculture","Saline Native","Saline Transition","Saline Monoculture"))

sample_data(datBac4c)$MarshClassVTransect<-factor(sample_data(datBac4c)$MarshClassVTransect,levels=c("Fresh Native","Fresh Transition","Fresh Phragmites","Brackish Native","Brackish Transition","Brackish Phragmites","Saline Native","Saline Transition","Saline Phragmites"))


tempphyB<-datBac4c%>%
  subset_samples(Transect!="Transition")%>%
#  subset_samples(Site=="Turtle Cove")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
anova(mynmdsB,by="margin",permutations = how(blocks=sample_data(tempphyB)$Site,nperm=9999))
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~MarshClassV+Transect))
anova(mynmdsB,by="margin",permutations = how(blocks=sample_data(tempphyB)$Site,nperm=9999))

#varianced explained
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsB
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsB
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~Transect*MarshClassV+Condition(Transect+MarshClassV)));mynmdsB

#it looks like the effect of phragmites varies by SITE rather than marsh class, the site*transect is much more significant (explains 10% of the variation)
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~Transect*Site+Condition(Transect+Site)));mynmdsB
anova(mynmdsB,by="terms",permutations = how(blocks=sample_data(tempphyB)$Site,nperm=9999))


#Plotting
mynmdsB <- ordinate(tempphyB, "CAP", "bray",formula=as.formula(~MarshClassVTransect))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/PCoAbacteria.pdf",width=4.8,height=3)
plot_ordination(tempphyB, mynmdsB, type="samples", color="MarshClassVTransect",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh Class/Invasion")+
  scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh Class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=MarshClass_Transect),level=.95)
dev.off()


##### Only on haplotype I #####
tempphyBI<-tempphyB%>%
  subset_samples(Site!="LUMCON 1",Site!="LUMCON 2")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmdsBI <- ordinate(tempphyBI, "CAP", "bray",formula=as.formula(~MarshClassV*Transect))
anova(mynmdsBI,by="terms",permutations = how(blocks=sample_data(tempphyBI)$Site,nperm=9999))
mynmdsBI <- ordinate(tempphyBI, "CAP", "bray",formula=as.formula(~MarshClassV+Transect))
anova(mynmdsBI,by="margin",permutations = how(blocks=sample_data(tempphyBI)$Site,nperm=9999))
mynmdsBI <- ordinate(tempphyBI, "CAP","bray",formula=as.formula(~MarshClassV+Condition(Transect)));mynmdsBI
mynmdsBI <- ordinate(tempphyBI, "CAP","bray",formula=as.formula(~Transect+Condition(MarshClassV)));mynmdsBI
mynmdsBI <- ordinate(tempphyBI, "CAP","bray",formula=as.formula(~MarshClassV*Transect+Condition(Transect+MarshClassV)));mynmdsBI



##### Bacteria - separate ordinations by marsh type #####
tempphyBF<-datBac4c%>%
  subset_samples(Plot!=23)%>% #take out outlier, not sure why it is different
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Fresh")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsBFresh <- ordinate(tempphyBF, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsBFresh,permutations = how(blocks=sample_data(tempphyBF)$Site,nperm=9999),by="margin")
mynmdsBFresh <- ordinate(tempphyBF, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Bac1<-plot_ordination(tempphyBF, mynmdsBFresh, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]



tempphyBB<-datBac4c%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Brackish")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsBBrackish <- ordinate(tempphyBB, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsBBrackish,permutations = how(blocks=sample_data(tempphyBB)$Site,nperm=9999))
mynmdsBBrackish <- ordinate(tempphyBB, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Bac2<-plot_ordination(tempphyBB, mynmdsBBrackish, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]



tempphyBS<-datBac4c%>%
  subset_samples(Transect!="Transition")%>%
  subset_samples(MarshClassV=="Saline")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)

mynmdsBSaline <- ordinate(tempphyBS, "CAP","bray",formula=as.formula(~Transect+Condition(Site)))#
anova(mynmdsBSaline,permutations = how(blocks=sample_data(tempphyBS)$Site,nperm=9999),by="margin")
#mynmdsBSaline <- ordinate(tempphyBS, "CAP","bray",formula=as.formula(~1+Condition(Site)))#
Bac3<-plot_ordination(tempphyBS, mynmdsBSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

#sppscores(mynmdsFFresh) <- data.frame(t(otu_table(datPlantsF)))
#plot(mynmdsFresh)
ind<-order(scores(mynmdsFFresh, display = "species")[,1])
scores(mynmdsFFresh, display = "species")[ind,1]


legendPlant<-plot_ordination(datPlantsS, mynmdsSaline, type=c("sites"), color="Transect",axes=c(1,2))+
  theme_classic()+
  geom_point(size = 2)+  
  scale_color_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  scale_fill_manual(values = c("#56ae6c", "#8960b3"),labels = c("Native", "Phragmites"),name = "Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Transect),level=.95)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/SeparateOrdinationsBaci.pdf",width=12.4,height=3)#4.8
legend <- cowplot::get_legend(legendPlant)
plot_grid(Bac1,Bac2,Bac3,legend,nrow = 1)
dev.off()





###### Percent explained of plants, fungi, bacteria##### 
percentexplained2<-percentexplained%>%
  gather(type,value,percentplant:sigbac)%>%
  mutate(taxon=rep(c("Plant","Fungi","Bacteria"),each=24))%>%
  mutate(type=rep(rep(c("percent","f","p"),each=8),3))%>%
  mutate(Site=rep(sites,9))%>%
  spread(type,value)
percentexplained2$Site<-factor(percentexplained2$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
percentexplained2$taxon<-factor(percentexplained2$taxon,levels=c("Plant","Fungi","Bacteria"))
percentexplained2<-percentexplained2%>%
  arrange(taxon,Site)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/percentexplained.pdf")
ggplot(percentexplained2,aes(x=Site,y=percent,group=taxon,fill=taxon))+
  labs(x = "",y="Phrag impact")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_bar(stat="identity",alpha=ifelse(percentexplained2$p<.05,1,.5))+
  facet_wrap(~taxon,nrow=3,scales="free")#
dev.off()

#one ordination (nmds) of bacteria data to see where pearl river groups. pearl river is right on top of big branch, next to turtle cove and then bayou sauvage. so it is grouping more stongly with veg
tempphy<-datBac4c%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect*Site))
mynmds <- ordinate(tempphy, "NMDS", "bray")
anova(mynmds,by="terms")
anova(mynmds,by="margin")
plot_ordination(tempphy, mynmds, type="samples", color="Site")+
  theme_classic()#+  theme(legend.position = "none")



#combinding ordination data together
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")
percentexplained3<-percentexplained
percentexplained3$Salinity<-factor(c("l","l","m","m","m","m","h","h"))
m1<-gls(percentplant~Salinity,data=percentexplained3)
anova(m1,type="marginal")
m1<-gls(percentfungi~Salinity,data=percentexplained3)
anova(m1,type="marginal")
m1<-gls(percentbac~Salinity,data=percentexplained3)
anova(m1,type="marginal")

percentexplained4<-percentexplained2%>%
  filter(taxon!="Plant")
percentexplained4$Salinity<-factor(rep(c("l","l","m","m","m","m","h","h")))

m1<-gls(percent~taxon*Salinity,data=percentexplained4)
anova(m1,type="marginal")

plot(percentexplained3$percentfungi,percentexplained3$percentbac)
#m1<-gls(percentfungi~percentbac,data=percentexplained3)
#anova(m1,type="marginal")
rcorr(percentexplained3$percentfungi,percentexplained3$percentbac) #same statistic as gls but I like it b/c it gives the correlation coefficient

plot(percentexplained3$percentplant,percentexplained3$percentfungi)
rcorr(percentexplained3$percentplant,percentexplained3$percentfungi)  

plot(percentexplained3$percentplant,percentexplained3$percentbac)
rcorr(percentexplained3$percentplant,percentexplained3$percentbac)  





#ordination code just for native and phrag plots
ordlistp<-list(NA)
percentexplained<-data.frame(percentplant=rep(NA,8),fplant=rep(NA,8),sigplant=rep(NA,8),percentfungi=rep(NA,8),ffungi=rep(NA,8),sigfungi=rep(NA,8),percentbac=rep(NA,8),fbac=rep(NA,8),sigbac=rep(NA,8))

for (i in 1:8){
  tempphy<-datPlants%>%
    subset_samples(Site==sites[i])%>%
    subset_samples(Transect!="Transition")%>%
    filter_taxa(function(x) sum(x>0) >0, prune=T)
  tempphy<-prune_samples(sample_sums(tempphy)>0,tempphy)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistp[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,1]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,2]<-tempanova$F[1]
  percentexplained[i,3]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistp[[1]],ordlistp[[2]],ordlistp[[3]],ordlistp[[4]],ordlistp[[5]],ordlistp[[6]],ordlistp[[7]],ordlistp[[8]],legend, nrow = 1)


#FUNGI
ordlistf<-list(NA)
for (i in 1:8){
  tempphy<-datITSS5c%>%
    subset_samples(Site==sites[i])%>%
    subset_samples(Transect!="Transition")%>%
    filter_taxa(function(x) sum(x>0) >2, prune=T)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistf[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,4]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,5]<-tempanova$F[1]
  percentexplained[i,6]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistf[[1]],ordlistf[[2]],ordlistf[[3]],ordlistf[[4]],ordlistf[[5]],ordlistf[[6]],ordlistf[[7]],ordlistf[[8]],legend, nrow = 1)

#BACTERIA
ordlistb<-list(NA)
for (i in 1:8){
  tempphy<-datBac4c%>%
    subset_samples(Plot!=168)%>% #take out crazy 2 million read plot
    subset_samples(Site==sites[i])%>%
    subset_samples(Transect!="Transition")%>%
    filter_taxa(function(x) sum(x>0) >2, prune=T)
  mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect))
  ordlistb[[i]]<-plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
    theme_classic()+
    theme(legend.position = "none")+
    stat_ellipse(geom = "polygon", type="t", alpha=0.4, aes(fill=Transect),level=.95)
  percentexplained[i,7]<-mynmds$CCA$tot.chi/mynmds$tot.chi
  tempanova<-anova(mynmds,by = "margin",permutations = how(nperm=10000))
  percentexplained[i,8]<-tempanova$F[1]
  percentexplained[i,9]<-tempanova$P[1]
}
legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
plot_grid(ordlistb[[1]],ordlistb[[2]],ordlistb[[3]],ordlistb[[4]],ordlistb[[5]],ordlistb[[6]],ordlistb[[7]],ordlistb[[8]],legend, nrow = 1)

percentexplained2<-percentexplained%>%
  gather(type,value,percentplant:sigbac)%>%
  mutate(taxon=rep(c("Plant","Fungi","Bacteria"),each=24))%>%
  mutate(type=rep(rep(c("percent","f","p"),each=8),3))%>%
  mutate(Site=rep(sites,9))%>%
  spread(type,value)
percentexplained2$Site<-factor(percentexplained2$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
percentexplained2$taxon<-factor(percentexplained2$taxon,levels=c("Plant","Fungi","Bacteria"))
percentexplained2<-percentexplained2%>%
  arrange(taxon,Site)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/percentexplainednativephrag.pdf")
ggplot(percentexplained2,aes(x=Site,y=percent,group=taxon,fill=taxon))+
  labs(x = "",y="Phrag impact")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_bar(stat="identity",alpha=ifelse(percentexplained2$p<.05,1,.5))+
  facet_wrap(~taxon,nrow=3,scales="free")#
dev.off()









##### Bray curtis / sorenson / jaccard in plot pairs (lines) #####

#Plants, not including phragmites
#abun<-dat17[,6:51]
#relabun<-dat17[,6:51]/rowSums(dat17[,6:51])

#Plants, including phragmites
abun<-dat17[,5:51]
relabun<-dat17[,5:51]/rowSums(dat17[,5:51])

#delete species that are not present in 2017
relabun2<-relabun[,which(colSums(relabun,na.rm=T)>0)]
#temp<-vegdist(relabun2,method="bray")
#as.matrix(temp)

output<-data.frame(Line=1:56,bray=NA,sor=NA,jac=NA)
for(i in 1:length(unique(dat17$Line))){
  indp<-which(dat17$Line==i&dat17$Transect=="Phragmites")
  indn<-which(dat17$Line==i&dat17$Transect=="Native")
  tempp<-relabun2[indp,]
  tempn<-relabun2[indn,]
  temppa<-abun[indp,]
  tempna<-abun[indn,]
  try(output$bray[i]<-vegdist(rbind(tempp,tempn),method="bray")) #dissimilarity
  output$sor[i]<-1-sorenson(temppa,tempna) #similarity in pres/abs
  try(output$jac[i]<-vegdist(rbind(temppa,tempna),method="jaccard",binary=T))
  output$Site[i]<-as.character(dat17$Site[indp])
}
output$Site<-factor(output$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
outputPlant<-output

#not including phrag should make sorenson's/jaccard denominator smaller and thus the similarity metric larger, and the disimilarity metric smaller

outputPlant2 <- outputPlant%>%
  group_by(Site)%>%
  summarise(mean=mean(bray,na.rm=T),se=std.error(bray,na.rm=T),sd=sd(bray,na.rm=T))
#outputPlant2 <- outputPlant%>%
#  group_by(Site)%>%
#  summarise(mean=mean(sor,na.rm=T),se=std.error(sor,na.rm=T),sd=sd(sor,na.rm=T))
outputPlant2 <- outputPlant%>%
  group_by(Site)%>%
  summarise(mean=mean(jac,na.rm=T),se=std.error(jac,na.rm=T),sd=sd(sor,na.rm=T))

#The fossile package returns similarity on all indexes, so need to do 1-returnvalue. it uses binary version of jaccard.
#If I don't include phragmites, bray does not make sense b/c it is measuring the difference in relative abundance of species in the plots (it doesn't make logical sense for plots with no species). Jaccard and Sorensons don't really make sense either but they can be computed, they will always be zero (or 1 if dissimilarity) if you compare to a plot wth no speices (because there are zero shared species)
#For including phrag analyses: we can calculate any of the metrics and they make sense. if you use bray, you don't get very much variation in plots from fontainbleau/big branch, bayou sauvage, LUM1, LUM2 b/c all of these had super dominant phrag in the monoculture plots and super dominant natives in native plots so they all look very very different from native plots. if you use jaccard or sorensons, then there is more variation since it is based on pres/abs not relative abundance. Sorensons and jaccard give very similar results

sorplant<-ggplot(outputPlant2,aes(x=Site,y=mean))+
  labs(x = "",y="Dissimilarity (Plant)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  #geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.8) 
sorplant

outputPlant3<-outputPlant%>%
  right_join(unique(dat17[,c("Site","MarshClassV","MarshClassS")]),"Site")

m1<-lme(jac~MarshClassV,random=~1|Site,data=outputPlant3)#varIdent(form =~ 1|Site), cant do this b/c some sites have no variance
anova(m1,type="marginal")
summary(glht(m1, linfct = mcp(MarshClassV = "Tukey")))

m1P<-as.data.frame(summary(emmeans(m1,~MarshClassV)))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/dissimilarityplant.pdf",width=2.2,height=2.2)
ggplot(m1P,aes(x=MarshClassV,y=emmean))+
  labs(x = "",y="Plant dissimilarity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5)
dev.off()




#FUNGI
sample_data(datITSS5c)

#A number of samples did not have reads
outputm<-data.frame(Line=1:56,bray=NA,sor=NA,jac=NA)
for(i in 1:56){
  temp<-sample_data(datITSS5c)
  indp<-which(temp$Line==i&temp$Transect=="Phragmites")
  indn<-which(temp$Line==i&temp$Transect=="Native")
  if(length(indp)+length(indn)==2){
    tempp<-otu_table(datITSS5c)[indp,]
    tempn<-otu_table(datITSS5c)[indn,]
    outputm$bray[i]<-vegdist(rbind(tempp,tempn),method="bray") 
    outputm$sor[i]<-1-sorenson(as.vector(tempp),as.vector(tempn)) 
    outputm$jac[i]<-vegdist(rbind(tempp,tempn),method="jaccard",binary=T)
    outputm$Site[i]<-as.character(temp$Site[indp])
  }else{
    outputm$bray[i]<-NA #dissimilarity
    outputm$sor[i]<-NA #similarity in pres/abs
    outputm$jac[i]<-NA #similarity in pres/abs
    outputm$Site[i]<-NA
  }  
}
outputm$Site<-factor(output$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
outputITS<-outputm

outputITS2 <- outputITS%>%
  group_by(Site)%>%
  summarise(mean=mean(bray,na.rm=T),se=std.error(bray,na.rm=T),sd=sd(bray,na.rm=T))
# outputITS2 <- outputITS%>%
#   group_by(Site)%>%
#   summarise(mean=mean(jac,na.rm=T),se=std.error(jac,na.rm=T),sd=sd(jac,na.rm=T))

sorITS<-ggplot(outputITS2,aes(x=Site,y=mean))+
  labs(x = "",y="Dissimilarity (Fungi)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  #geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.8) 
sorITS

outputITS3<-outputITS%>%
  right_join(unique(dat17[,c("Site","MarshClassV","MarshClassS")]),"Site")

m1<-lme(bray~MarshClassV,random=~1|Site,varIdent(form =~ 1|Site),data=outputITS3,na.action=na.omit)
anova(m1,type="marginal")

m1F<-as.data.frame(summary(emmeans(m1,~MarshClassV)))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/dissimilarityfungi.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=MarshClassV,y=emmean))+
  labs(x = "",y="Fungi dissimilarity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5)
dev.off()



#Checking, how does within native and within phrag diversity compare to between phrag-native. Within native and within phrag beta diveristy varies substiantialy across sites, thus the dissimilarity plots above are strange b/c you are not accounting for different overal betadiversity within a site 
# temp<-sample_data(datITSS5c)
# indp<-which(temp$Site=="Barataria"&temp$Transect=="Phragmites")
# indn<-which(temp$Site=="Barataria"&temp$Transect=="Native")
# tempp<-otu_table(datITSS5c)[indp,]
# tempn<-otu_table(datITSS5c)[indn,]
# mean(vegdist(tempp,method="bray")) 
# mean(vegdist(tempn,method="bray")) 




#Bacteria
datBac4c

#A number of samples did not have reads
outputb<-data.frame(Line=1:56,bray=NA,sor=NA,jac=NA)
#remove plots 168 with 2 million reads
datBac4cno168<-datBac4c%>%
  subset_samples(Plot!=168) #take out crazy 2 million read plot

for(i in 1:56){
  temp<-sample_data(datBac4cno168)
  indp<-which(temp$Line==i&temp$Transect=="Phragmites")
  indn<-which(temp$Line==i&temp$Transect=="Native")
  if(length(indp)+length(indn)==2){
    tempp<-t(otu_table(datBac4cno168))[indp,]
    tempn<-t(otu_table(datBac4cno168))[indn,]
    outputb$bray[i]<-vegdist(rbind(tempp,tempn),method="bray") 
    outputb$sor[i]<-1-sorenson(as.vector(tempp),as.vector(tempn)) 
    outputb$jac[i]<-vegdist(rbind(tempp,tempn),method="jaccard",binary=T)
    outputb$Site[i]<-as.character(temp$Site[indp])
  }else{
    outputb$bray[i]<-NA #dissimilarity
    outputb$sor[i]<-NA #similarity in pres/abs
    outputb$jac[i]<-NA #similarity in pres/abs
    outputb$Site[i]<-NA
  }  
}
outputb$Site<-factor(output$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
outputBac<-outputb

outputBac2 <- outputBac%>%
  group_by(Site)%>%
  summarise(mean=mean(bray,na.rm=T),se=std.error(bray,na.rm=T),sd=sd(bray,na.rm=T))
# outputBac2 <- outputBac%>%
#   group_by(Site)%>%
#   summarise(mean=mean(jac,na.rm=T),se=std.error(jac,na.rm=T),sd=sd(jac,na.rm=T))

sorbac<-ggplot(outputBac2,aes(x=Site,y=mean))+
  labs(x = "",y="Dissimilarity (Bacteria)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.8) 
sorbac

outputBac3<-outputBac%>%
  right_join(unique(dat17[,c("Site","MarshClassV","MarshClassS")]),"Site")

m1<-lme(bray~MarshClassV,random=~1|Site,varIdent(form =~ 1|Site),data=outputBac3,na.action=na.omit)
#m1<-gls(bray~MarshClassV,data=outputBac3,na.action=na.omit)
anova(m1,type="marginal")
summary(glht(m1, linfct = mcp(MarshClassV = "Tukey")))


m1B<-as.data.frame(summary(emmeans(m1,~MarshClassV)))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/dissimilarityfungi.pdf",width=2.2,height=2.2)
ggplot(m1B,aes(x=MarshClassV,y=emmean))+
  labs(x = "",y="Bacteria dissimilarity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5)
dev.off()


#calcualting disimlarity for both phrag and transition plots and rbinding them
outputBactot<-rbind(outputBac,outputb)
outputBactot$Transect<-rep(c("Phragmites","Transition"),each=56)
outputBactot3<-outputBactot%>%
  left_join(unique(dat17[,c("Site","MarshClassV","MarshClassS")]),"Site")
outputBactot3$MarshClassVTransect<-interaction(outputBactot3$MarshClassV,outputBactot3$Transect)
m1<-lme(bray~MarshClassVTransect,random=~1|Site,data=outputBactot3,na.action=na.omit)
#m1<-gls(bray~MarshClassV,data=outputBac3,na.action=na.omit)
anova(m1,type="marginal")
summary(glht(m1, linfct = mcp(MarshClassVTransect = "Tukey")))

m1B<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))
ggplot(m1B,aes(x=Transect,y=emmean,group=MarshClassV))+
  labs(x = "",y="Bacteria dissimilarity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5)+
  facet_wrap(~MarshClassV)
#




#plot of plants, fungi, bacteria
#legend <- cowplot::get_legend(plot_ordination(tempphy, mynmds, type="samples", color="Transect"))
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/braydissimilarityincphrag.pdf")#,width=6,height=4
plot_grid(sorplant,sorITS,sorbac,nrow = 3)
dev.off()




##### Coordinated Response Analysis - Composition #####

outputPlant2
outputITS2
outputBac2


outputPlant2$sd2<-outputPlant2$sd
outputPlant2$sd2[c(4,8)]<-.01

m1<-deming(outputITS2$mean~outputPlant2$mean,xstd=outputPlant2$sd2,ystd = outputITS2$sd)#,cv=T or F gives same exact results
m1

m1<-deming(outputBac2$mean~outputPlant2$mean,xstd=outputPlant2$sd2,ystd = outputBac2$sd)#,cv=T or F gives same exact results
m1

m1<-deming(outputBac2$mean~outputITS2$mean,xstd=outputITS2$sd,ystd = outputBac2$sd,conf=.95)#,cv=T or F gives same exact results
m1


#Only haplotype I
outputPlant2I<-outputPlant2[1:6,]
outputITS2I<-outputITS2[1:6,]
outputBac2I<-outputBac2[1:6,]

m1<-deming(outputITS2I$mean~outputPlant2I$mean,xstd=outputPlant2I$sd2,ystd = outputITS2I$sd)#,cv=T or F gives same exact results
m1

m1<-deming(outputBac2I$mean~outputPlant2I$mean,xstd=outputPlant2I$sd2,ystd = outputBac2I$sd)#,cv=T or F gives same exact results
m1

m1<-deming(outputBac2I$mean~outputITS2I$mean,xstd=outputITS2I$sd,ystd = outputBac2I$sd,conf=.95)#,cv=T or F gives same exact results
m1



#plotting 

outputAll<-data.frame(outputPlant,distfungi=outputITS$bray,distbac=outputBac$bray)
colnames(outputAll)[c(2,3)]<-c("distplant","sorplant")

outputAll2<-outputAll%>%
  dplyr::select(Site,jac,distfungi,distbac)%>%
  group_by(Site)%>%
  summarise_all(list(~mean(.,na.rm=T),~sd(.,na.rm=T),~std.error(.,na.rm=T)))
as.data.frame(outputAll2)


pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/compdifplantfungi.pdf",width=2,height=1.9)#2.2
ggplot(outputAll2,aes(x=jac_mean,y=distfungi_mean))+
  labs(x = "Composition difference plant",y="Composition difference fungi")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.4))+
  geom_point(size=1.6)+
  geom_errorbar(aes(ymax = distfungi_mean+distfungi_std.error, ymin=distfungi_mean-distfungi_std.error),width=0,size=.4)+#
  geom_errorbarh(aes(xmax=jac_mean+jac_std.error,xmin=jac_mean-jac_std.error),height=0,size=.4)+#
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.2466269,slope=0.6213904,size=.4)
dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/compdifplantbac.pdf",width=2,height=1.9)
ggplot(outputAll2,aes(x=jac_mean,y=distbac_mean))+
  labs(x = "Composition difference plant",y="Composition difference bacteria")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.4))+
  geom_point(size=1.6)+
  geom_errorbar(aes(ymax = distbac_mean+distbac_std.error, ymin=distbac_mean-distbac_std.error),width=0,size=.4)+#
  geom_errorbarh(aes(xmax=jac_mean+jac_std.error,xmin=jac_mean-jac_std.error),height=0,size=.4)+#
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.2058173,slope=0.6072133,size=.4)
dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/compdiffungibac.pdf",width=2,height=1.9)
ggplot(outputAll2,aes(x=distfungi_mean,y=distbac_mean))+
  labs(x = "Composition difference fungi",y="Composition difference bacteria")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.4))+
  geom_point(size=1.6)+
  geom_errorbar(aes(ymax = distbac_mean+distbac_std.error, ymin=distbac_mean-distbac_std.error),width=0,size=.4)+#
  geom_errorbarh(aes(xmax=distfungi_mean+distfungi_std.error,xmin=distfungi_mean-distfungi_std.error),height=0,size=.4)+#
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.01262485,slope=0.91522253,size=.4)
dev.off()


#old analysis of all points

plot(outputPlant$sor,outputITS$dist)
abline(lm(outputITS$dist~outputPlant$sor))
#to plot this corectly I probably have to use predict (level=0)
m1<-lme(distfungi~sor,random=~1|Site,data=outputAll,na.action=na.omit)
anova(m1)
plot(outputPlant2$mean,outputITS2$mean)
abline(m1,col=3)
abline(lm(outputITS2$mean~outputPlant2$mean),col=2)
plot(outputPlant2$mean,outputBac2$mean)
plot(outputITS2$mean,outputBac2$mean)

anova(lm(outputPlant2$mean~outputITS2$mean))
rcorr(outputPlant2$mean,outputITS2$mean)
rcorr(outputPlant2$mean,outputBac2$mean)
rcorr(outputITS2$mean,outputBac2$mean)






###### FUNguild #####

##### Pathotroph #####
#m1<-lme(Pathotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(Pathotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Pathotroph~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Pathogen percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()


##### Symbiotroph #####
#m1<-lme(Symbiotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(Symbiotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)#,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T)
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Symbiotroph~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Pathogen percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Saprotroph #####
m1<-lme(Saprotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m2<-lme(Saprotroph~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)#,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T)
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Saprotroph~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Saprotroph percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Plant Pathogen #####

#looking at outliers in saline native and trans plots
plot(sqrt(dat17$Plant.Pathogen)~dat17$MarshClassV.Transect)
boxplot(dat17$Plant.Pathogen)

dat17$Plant.Pathogen.nooutliers<-dat17$Plant.Pathogen
dat17$Plant.Pathogen.nooutliers[which(dat17$Plant.Pathogen.nooutliers>10)]<-NA
boxplot(dat17$Plant.Pathogen.nooutliers)

#m1<-lme(sqrt(Plant.Pathogen.nooutliers)~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(log(Plant.Pathogen.nooutliers+1)~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

#make figure with non-transformed data??
#m1a<-lme(Plant.Pathogen.nooutliers~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)#,control = lmeControl(msMaxIter = 10000,tolerance=1e-06
#m1Fa<-as.data.frame(summary(emmeans(m1a,~Transect|MarshClassV)))

m2<-lme(log(Plant.Pathogen.nooutliers+1)~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#doing glht with interactions the "correct" way with m1, it gives the same exact result as above
tmp <- expand.grid(Transect = factor(levels(dat17$Transect),levels=c("Native","Transition","Phragmites")),MarshClassV = factor(levels(dat17$MarshClassV),levels=c("Fresh","Brackish","Saline")))
X <- model.matrix(~ MarshClassV * Transect, data = tmp)
glht(m1, linfct = X)
predict(m1, newdata = tmp,level=0)
Tukey <- contrMat(table(dat17$Transect), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K1) <- paste(levels(dat17$MarshClassV)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K2) <- paste(levels(dat17$MarshClassV)[2], rownames(K2), sep = ":")
K3 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)),matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
rownames(K3) <- paste(levels(dat17$MarshClassV)[3], rownames(K3), sep = ":")
K <- rbind(K1, K2,K3)
colnames(K) <- c(colnames(Tukey), colnames(Tukey), colnames(Tukey))
summary(glht(m1, linfct = K %*% X))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Only on fresh and brackish, like carolyn's paper. results are smilar but not exactly the same, likely due to her 50% (not 70%) confidence in taxonomy assignment
m1<-lme(Plant.Pathogen.nooutliers~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=subset(dat17,dat17$MarshClassV!="Saline"),na.action=na.omit)
m2<-lme(Plant.Pathogen.nooutliers~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=subset(dat17,dat17$MarshClassV!="Saline"),na.action=na.omit)
anova(m1,m2) #het var sig
anova(m2,type="marginal")
summary(glht(m2,linfct=mcp(Transect="Tukey")))
m1F<-as.data.frame(summary(emmeans(m2,~Transect)))

#fig
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant pathogen percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()

#Make figure by back-transforming
m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))
m1Ftrans<-m1F[,1:2]
m1Ftrans$emmean<-exp(m1F$emmean)-1
m1Ftrans$lower<-exp(m1F$emmean-m1F$SE)-1
m1Ftrans$upper<-exp(m1F$emmean+m1F$SE)-1
#log(data+1)=newdata
#data=exp(newdata)-1

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/pathogenabundance.pdf",width=2.2,height=2.2)
ggplot(m1Ftrans,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant pathogen percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = upper, ymin=lower),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Arbuscular.Mycorrhizal #####

#the models with fixed effects MarshClassV*Transect vs MarshClassV.Transect are only identical if fit with ML

hist(dat17$Arbuscular.Mycorrhizal)
plot(log(dat17$Arbuscular.Mycorrhizal+1)~dat17$MarshClassV.Transect)

#m1<-lme(log(Arbuscular.Mycorrhizal+1)~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(log(Arbuscular.Mycorrhizal+1)~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))
plot(dat17$Transect[is.na(dat17$Arbuscular.Mycorrhizal)==F],resid(m1,type="normalized"))
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Arbuscular.Mycorrhizal~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

##Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="AMF percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()

#Make figure by back-transforming
m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))
m1Ftrans<-m1F[,1:2]
m1Ftrans$emmean<-exp(m1F$emmean)-1
m1Ftrans$lower<-exp(m1F$emmean-m1F$SE)-1
m1Ftrans$upper<-exp(m1F$emmean+m1F$SE)-1
#log(data+1)=newdata
#data=exp(newdata)-1

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/AMFabundance.pdf",width=2.2,height=2.2)
ggplot(m1Ftrans,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="AMF percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = upper, ymin=lower),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Ectomycorrhizal - too low abundance to be that useful, not zeros but just very low abundance #####
plot(dat17$Ectomycorrhizal~dat17$MarshClassV.Transect)
boxplot(dat17$Ectomycorrhizal)

dat17$Ectomycorrhizal.nooutliers<-dat17$Ectomycorrhizal
dat17$Ectomycorrhizal.nooutliers[which(dat17$Ectomycorrhizal.nooutliers>1.5)]<-NA
boxplot(dat17$Arbuscular.Mycorrhizal)

#m1<-lme(Ectomycorrhizal~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(Ectomycorrhizal.nooutliers~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Ectomycorrhizal~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

##Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Ecto percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()




##### PlantMutualist = Ectomycorrhizal+Arbuscular.Mycorrhizal+Endophyte #####
plot(dat17$PlantMutualist~dat17$MarshClassV.Transect)
boxplot(dat17$PlantMutualist)

#m1<-lme(PlantMutualist~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(PlantMutualist~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(PlantMutualist~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

##Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant mutualist percent abundance") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Endophyte - too many zeros #####
##### Plant.Pathogenanywhere #####

plot(dat17$Plant.Pathogenanywhere~dat17$MarshClassV.Transect)
hist(dat17$Plant.Pathogenanywhere)

#m1<-lme(Plant.Pathogenanywhere~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(Plant.Pathogenanywhere~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Plant.Pathogenanywhere~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Only on fresh and brackish, like carolyn's paper. 
carodat<-dat17%>%
  filter(MarshClassV!="Saline")%>%
  filter(Transect!="Transition")%>%
  filter(Plot!=88)
plot(carodat$Plant.Pathogenanywhere~carodat$MarshClassV.Transect)
carodat$Plant.Pathogenanywherenoout<-carodat$Plant.Pathogenanywhere
carodat$Plant.Pathogenanywherenoout[which(carodat$Plant.Pathogenanywherenoout>15)]<-NA

m2<-lme(Plant.Pathogenanywherenoout~Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|Transect),data=carodat,na.action=na.omit)
anova(m2,type="marginal")
summary(glht(m2,linfct=mcp(Transect="Tukey")))
m1F<-as.data.frame(summary(emmeans(m2,~Transect)))
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=Transect))+
  labs(x = "",y="Plant pathogen anywhere") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))


##Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant pathogen anywhere") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()
