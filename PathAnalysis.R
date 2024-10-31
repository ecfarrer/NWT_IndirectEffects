#Path analysis


##### NMDS outputs for path analysis #####

datITSS5cE
dat16SS5cE

###### Wet meadow 16S ######

tempphyEWM16S<-dat16SS5cE%>%
  subset_samples(SampleType=="roots")%>%
  subset_samples(CommunityType=="WM")%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x))  #note that filtering >2 can take out as many has 3800 reads (from 7198)
#sample_data(tempphyEWM)$CommunityPlotType<-paste(sample_data(tempphyEWM)$CommunityType,sample_data(tempphyEWM)$PlotType,sep="")

rowSums(otu_table(tempphyEWM16S))

mynmdsE <- ordinate(tempphyEWM16S, "NMDS", "bray")
#mynmdsE <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~PlotType))
#anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWM16S, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyEWM16S, mynmdsE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)

mynmdsE <- ordinate(tempphyEWM16S, "CAP", "bray",formula=as.formula(~1))
scores(mynmdsE)$sites
sample_data(tempphyEWM16S)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWM16S)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEWM16S, "NMDS", "bray")
sample_data(tempphyEWM16S)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWM16S)$NMDS2<-scores(mynmdsE)$sites[,2]

ind<-which(colnames(sample_data(tempphyEWM16S))=="MDS1")
colnames(sample_data(tempphyEWM16S))[ind]<-"bacteriaMDS1"
ind<-which(colnames(sample_data(tempphyEWM16S))=="MDS2")
colnames(sample_data(tempphyEWM16S))[ind]<-"bacteriaMDS2"
ind<-which(colnames(sample_data(tempphyEWM16S))=="NMDS1")
colnames(sample_data(tempphyEWM16S))[ind]<-"bacteriaNMDS1"
ind<-which(colnames(sample_data(tempphyEWM16S))=="NMDS2")
colnames(sample_data(tempphyEWM16S))[ind]<-"bacteriaNMDS2"

sample_data(tempphyEWM16S)<-sample_data(tempphyEWM16S)[,cbind("PlotID","Treatment","bacteriaMDS1","bacteriaMDS2","bacteriaNMDS1","bacteriaNMDS2")]


boxplot(sample_data(tempphyEWM16S)$NMDS2~sample_data(tempphyEWM16S)$Treatment)

sample_data(tempphyEWM16S)



###### Wet meadow ITS ######

tempphyEWMITS<-datITSS5cE%>%
  subset_samples(SampleType=="roots")%>%
  subset_samples(CommunityType=="WM")%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x))  #note that filtering >2 can take out as many has 3800 reads (from 7198)
#sample_data(tempphyEWM)$CommunityPlotType<-paste(sample_data(tempphyEWM)$CommunityType,sample_data(tempphyEWM)$PlotType,sep="")

rowSums(otu_table(tempphyEWMITS))

mynmdsE <- ordinate(tempphyEWMITS, "NMDS", "bray")
#mynmdsE <- ordinate(tempphyEWM, "CAP", "bray",formula=as.formula(~PlotType))
#anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWMITS, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyEWMITS, mynmdsE, type="samples", color="Treatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Treatment),level=.95)

mynmdsE <- ordinate(tempphyEWMITS, "CAP", "bray",formula=as.formula(~1))
scores(mynmdsE)$sites
sample_data(tempphyEWMITS)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMITS)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEWMITS, "NMDS", "bray")
sample_data(tempphyEWMITS)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMITS)$NMDS2<-scores(mynmdsE)$sites[,2]

ind<-which(colnames(sample_data(tempphyEWMITS))=="MDS1")
colnames(sample_data(tempphyEWMITS))[ind]<-"fungiMDS1"
ind<-which(colnames(sample_data(tempphyEWMITS))=="MDS2")
colnames(sample_data(tempphyEWMITS))[ind]<-"fungiMDS2"
ind<-which(colnames(sample_data(tempphyEWMITS))=="NMDS1")
colnames(sample_data(tempphyEWMITS))[ind]<-"fungiNMDS1"
ind<-which(colnames(sample_data(tempphyEWMITS))=="NMDS2")
colnames(sample_data(tempphyEWMITS))[ind]<-"fungiNMDS2"

sample_data(tempphyEWMITS)<-sample_data(tempphyEWMITS)[,cbind("PlotID","Treatment","fungiMDS1","fungiMDS2","fungiNMDS1","fungiNMDS2")]

boxplot(sample_data(tempphyEWMITS)$MDS1~sample_data(tempphyEWMITS)$Treatment)

sample_data(tempphyEWMITS)



###### Wet and moist meadow 16S #####

tempphyEWMMM16S<-dat16SS5cE%>%
  subset_samples(SampleType=="roots")%>%
  subset_samples(CommunityType%in%c("WM","MM"))%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x)) 

rowSums(otu_table(tempphyEWMMM16S))

mynmdsE <- ordinate(tempphyEWMMM16S, "NMDS", "bray")
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~SiteTreatment))
anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyEWMMM16S, mynmdsE, type="samples", color="SiteTreatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=SiteTreatment),level=.95)

mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~1))
scores(mynmdsE)$sites
sample_data(tempphyEWMMM16S)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMM16S)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEWMMM16S, "NMDS", "bray")
sample_data(tempphyEWMMM16S)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMM16S)$NMDS2<-scores(mynmdsE)$sites[,2]

ind<-which(colnames(sample_data(tempphyEWMMM16S))=="MDS1")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS1"
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="MDS2")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS2"
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="NMDS1")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaNMDS1"
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="NMDS2")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaNMDS2"

sample_data(tempphyEWMMM16S)<-sample_data(tempphyEWMMM16S)[,cbind("PlotID","Treatment","bacteriaMDS1","bacteriaMDS2","bacteriaNMDS1","bacteriaNMDS2")]


boxplot(sample_data(tempphyEWMMM16S)$NMDS2~sample_data(tempphyEWMMM16S)$Treatment)

sample_data(tempphyEWMMM16S)






##### Path analysis #####

###### Wet meadow ######

datEWM<-datE%>%
  filter(CommunityType=="WM")
tempphyEWMITS
tempphyEWM16S
# temp<-sample_data(tempphyEWM)
# temp$SampleID<-gsub("_", "-", temp$PlotID)
# temp$SampleID<-sub("(-.*?)-", "\\1", temp$SampleID)

datEWMb<-datEWM%>%
  full_join(sample_data(tempphyEWM16S))%>%
  full_join(sample_data(tempphyEWMITS),by=c("PlotID","Treatment"))%>%
  filter(!is.na(bacteriaMDS1))%>%
  filter(!is.na(fungiMDS1))

#do i need to pre-standardize everything? yes it looks like particularly due to the Proj binary variable which is not standardized or centered (?) using the std.ov=T
datEWMc<-datEWMb%>%
  dplyr::select(Treatment, bacteriaNMDS1, bacteriaNMDS2, bacteriaMDS1, bacteriaMDS2,fungiNMDS1, fungiNMDS2, fungiMDS1, fungiMDS2, Biomassg,LeafNumber,LeafLengthmm, FlowersperRosette, Ecorrectedredo, Acorrectedredo,gsw)%>%
  dplyr::rename(Transpiration=Ecorrectedredo,Assimilation=Acorrectedredo)%>%
  mutate(Treatment=ifelse(Treatment=="Control",0,1))

datEWMcstand<-data.frame(apply(datEWMc,2,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))


#Using NMDS, I like NMDS b/c the points are not as far apart (like outliers) on the ordination axis (which makes some outlier points have too much leverage). I took out A corrected at the end b/c I realized that only growth/flowers was originally going to be in the path diagram
mod.wm = '
bacteriaNMDS1 ~ Treatment
bacteriaNMDS2 ~ Treatment
fungiNMDS1 ~ Treatment
fungiNMDS2 ~ Treatment
Biomassg ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
LeafLengthmm ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
LeafNumber ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
FlowersperRosette ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
Assimilation ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
Transpiration ~ Treatment + bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2
'

mod.wm = '
#bacteriaNMDS1 ~ Treatment
#bacteriaNMDS2 ~ Treatment
fungiNMDS1 ~ Treatment
fungiNMDS2 ~ Treatment
Biomassg ~  bacteriaNMDS1 + fungiNMDS1 + fungiNMDS2#Treatment ++ bacteriaNMDS2 
LeafLengthmm ~  bacteriaNMDS1  + fungiNMDS1 + fungiNMDS2#+ bacteriaNMDS2 Treatment +
LeafNumber ~ Treatment + bacteriaNMDS1  #+ fungiNMDS1 + fungiNMDS2 + bacteriaNMDS2
FlowersperRosette ~ Treatment+ fungiNMDS2   #+ bacteriaNMDS1 + bacteriaNMDS2 + fungiNMDS1
Assimilation ~ Treatment #+ fungiNMDS2+ bacteriaNMDS2  + fungiNMDS1 + bacteriaNMDS1
Transpiration ~ Treatment + bacteriaNMDS2 + fungiNMDS1 + fungiNMDS2 #+ bacteriaNMDS1#

Biomassg~~0*LeafLengthmm+0*LeafNumber+0*FlowersperRosette+0*Assimilation+0*Transpiration
LeafLengthmm~~0*FlowersperRosette+0*Assimilation+0*Transpiration
LeafNumber~~0*FlowersperRosette+0*Assimilation
FlowersperRosette~~0*Assimilation+0*Transpiration
'

mod.wm = sem(model = mod.wm, data = datEWMcstand) #datEWMb
summary(mod.wm, fit.measures=TRUE,rsquare=T)
AIC(mod.wm)

#Direct: 
0.508+1.210+0.588
#Indirect: 
-0.736*0.966+-0.736*1.023

boxplot(datEWMcstand$Transpiration~as.factor(datEWMcstand$Treatment))

ggplot(datEWMcstand, aes(x=LeafNumber, y=Transpiration, color=Treatment))+
  theme_classic()+
  geom_point()+
  geom_smooth(method="lm",se=F,color="black")

summary(lm(Transpiration~Treatment+bacteriaNMDS2+fungiNMDS1+fungiNMDS2,data=datEWMcstand))



###### Wet and moist meadow together ######











###### Dry meadow ######

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

