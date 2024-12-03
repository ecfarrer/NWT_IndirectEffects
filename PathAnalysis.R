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
#The 2 sites don't overlap in their composition. this will make the path analysis weird if i just use the straight NMDS b/c it will attribute site differences to microbes. things I can do: 1. put a site main effect in the path analysis and condition on site for the microbes (only can use MDS). 2. condition on site and don't put a site main effect (check regressions in path analysis). 3. condition on site and constrain by treatment (this might be good b/c the unconstrained by treatment ordinations never sort out by treatment so I will not get that arrow between treatment and microbes even though if i test just that relationship in an ordination it is significant)
#I tried conditioning on site, with no main effect of site and with the 2 bacterial axes the 2nd one was affected by climate change treatment (even though it didn't look like it from the ordination.). But now my concern is that if you get a strong effect of say bacteria axis 1 on leaf length and treatment does not affect bacteria axis 1, what does that really mean. It could be that bacterial axis 1 and leaf length are both responding to a nitrogen gradient or some other abiotic variable and that is not really a microbial affect. Because I don't have other abiotic measurements in the model I can't test this. I think this is a little sketchy.

tempphyEWMMM16S<-dat16SS5cE%>%
  subset_samples(SampleType=="roots")%>%
  subset_samples(CommunityType%in%c("WM","MM"))%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x)) 

rowSums(otu_table(tempphyEWMMM16S))

mynmdsE <- ordinate(tempphyEWMMM16S, "NMDS", "bray")
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~1+Condition(Site)))
anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyEWMMM16S, mynmdsE, type="samples", color="SiteTreatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=SiteTreatment),level=.95)

#unconstrained, conditioned on site
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~1+Condition(Site)))
scores(mynmdsE)$sites
sample_data(tempphyEWMMM16S)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMM16S)$MDS2<-scores(mynmdsE)$sites[,2]
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="MDS1")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS1uncsite"
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="MDS2")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS2uncsite"

#constrained by treatment, conditioned on site
mynmdsE <- ordinate(tempphyEWMMM16S, "CAP", "bray",formula=as.formula(~Treatment+Condition(Site)))
sample_data(tempphyEWMMM16S)$CAP1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMM16S)$MDS1<-scores(mynmdsE)$sites[,2]
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="CAP1")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS1tresite"
ind<-which(colnames(sample_data(tempphyEWMMM16S))=="MDS1")
colnames(sample_data(tempphyEWMMM16S))[ind]<-"bacteriaMDS2tresite"

sample_data(tempphyEWMMM16S)

sample_data(tempphyEWMMM16S)<-sample_data(tempphyEWMMM16S)[,cbind("PlotID","Treatment","bacteriaMDS1uncsite","bacteriaMDS2uncsite","bacteriaMDS1tresite","bacteriaMDS2tresite")]


boxplot(sample_data(tempphyEWMMM16S)$NMDS2~sample_data(tempphyEWMMM16S)$Treatment)

sample_data(tempphyEWMMM16S)



###### Wet and moist meadow ITS #####

tempphyEWMMMITS<-datITSS5cE%>%
  subset_samples(SampleType=="roots")%>%
  subset_samples(CommunityType%in%c("WM","MM"))%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)%>%
  transform_sample_counts(function(x) x/sum(x)) 

rowSums(otu_table(tempphyEWMMMITS))

mynmdsE <- ordinate(tempphyEWMMMITS, "NMDS", "bray")
mynmdsE <- ordinate(tempphyEWMMMITS, "CAP", "bray",formula=as.formula(~1+Condition(Site)))
mynmdsE <- ordinate(tempphyEWMMMITS, "CAP", "bray",formula=as.formula(~Treatment+Condition(Site)))
anova(mynmdsE)
mynmdsE <- ordinate(tempphyEWMMMITS, "CAP", "bray",formula=as.formula(~1))
plot_ordination(tempphyEWMMMITS, mynmdsE, type="samples", color="SiteTreatment",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=SiteTreatment),level=.95)

mynmdsE <- ordinate(tempphyEWMMMITS, "CAP", "bray",formula=as.formula(~1))
scores(mynmdsE)$sites
sample_data(tempphyEWMMMITS)$MDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMMITS)$MDS2<-scores(mynmdsE)$sites[,2]

mynmdsE <- ordinate(tempphyEWMMMITS, "NMDS", "bray")
sample_data(tempphyEWMMMITS)$NMDS1<-scores(mynmdsE)$sites[,1]
sample_data(tempphyEWMMMITS)$NMDS2<-scores(mynmdsE)$sites[,2]

ind<-which(colnames(sample_data(tempphyEWMMMITS))=="MDS1")
colnames(sample_data(tempphyEWMMMITS))[ind]<-"bacteriaMDS1"
ind<-which(colnames(sample_data(tempphyEWMMMITS))=="MDS2")
colnames(sample_data(tempphyEWMMMITS))[ind]<-"bacteriaMDS2"
ind<-which(colnames(sample_data(tempphyEWMMMITS))=="NMDS1")
colnames(sample_data(tempphyEWMMMITS))[ind]<-"bacteriaNMDS1"
ind<-which(colnames(sample_data(tempphyEWMMMITS))=="NMDS2")
colnames(sample_data(tempphyEWMMMITS))[ind]<-"bacteriaNMDS2"

sample_data(tempphyEWMMM16S)<-sample_data(tempphyEWMMM16S)[,cbind("PlotID","Treatment","bacteriaMDS1","bacteriaMDS2","bacteriaNMDS1","bacteriaNMDS2")]


boxplot(sample_data(tempphyEWMMM16S)$NMDS2~sample_data(tempphyEWMMM16S)$Treatment)

sample_data(tempphyEWMMM16S)





##### Path analysis #####

###### Wet meadow ######
#I'm worried about sample size. things I can try: 1. retain NAs in dataset, 2. combine WM and MM, 3. only use axis 1 for microbes to limit the number of variables

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


#[from annual report: Using NMDS, I like NMDS b/c the points are not as far apart (like outliers) on the ordination axis (which makes some outlier points have too much leverage). I took out A corrected at the end b/c I realized that only growth/flowers was originally going to be in the path diagram]
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

datEWMMM<-datE%>%
  filter(CommunityType%in%c("MM","WM"))
tempphyEWMMM16S
tempphyEWMMMITS
# temp<-sample_data(tempphyEWM)
# temp$SampleID<-gsub("_", "-", temp$PlotID)
# temp$SampleID<-sub("(-.*?)-", "\\1", temp$SampleID)

datEWMMMb<-datEWMMM%>%
  full_join(sample_data(tempphyEWMMM16S))%>%
  #full_join(sample_data(tempphyEWMITS))%>%
  filter(!is.na(bacteriaMDS1uncsite))%>%
  filter(!is.na(Ecorrectedredo))%>%
  filter(!is.na(FlowersperRosette))
  #filter(!is.na(fungiMDS1))

#do i need to pre-standardize everything? yes it looks like particularly due to the Proj binary variable which is not standardized or centered (?) using the std.ov=T
datEWMMMc<-datEWMMMb%>%
  dplyr::select(Treatment, bacteriaMDS1uncsite, bacteriaMDS2uncsite, bacteriaMDS1tresite, bacteriaMDS2tresite, Biomassg,LeafNumber,LeafLengthmm, FlowersperRosette, Ecorrectedredo, Acorrectedredo,gsw)%>%
  dplyr::rename(Transpiration=Ecorrectedredo,Assimilation=Acorrectedredo)%>%
  mutate(Treatment=ifelse(Treatment=="Control",0,1))

datEWMMMcstand<-data.frame(apply(datEWMMMc,2,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))
datEWMMMcstand$Site<-datEWMMMb$Site

#[from annual report: Using NMDS, I like NMDS b/c the points are not as far apart (like outliers) on the ordination axis (which makes some outlier points have too much leverage). I took out A corrected at the end b/c I realized that only growth/flowers was originally going to be in the path diagram]
mod.wmmm = '
bacteriaMDS1uncsite ~ Treatment
bacteriaMDS2uncsite ~ Treatment
Biomassg ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
LeafLengthmm ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
LeafNumber ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
FlowersperRosette ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
Assimilation ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
Transpiration ~ Treatment + bacteriaMDS1uncsite + bacteriaMDS2uncsite
'

#unconstrained, conditioned on site
mod.wmmm = '
#bacteriaMDS1uncsite ~ Treatment
bacteriaMDS2uncsite ~ Treatment
#Biomassg ~   bacteriaMDS2uncsite #Treatment +bacteriaMDS1uncsite +
LeafLengthmm ~ Treatment + bacteriaMDS2uncsite + bacteriaMDS1uncsite# 
#LeafNumber ~ bacteriaMDS1uncsite  #+ bacteriaMDS2uncsite Treatment + 
FlowersperRosette ~ bacteriaMDS1uncsite #+ bacteriaMDS2uncsite #Treatment + 
#Assimilation ~  bacteriaMDS1uncsite #+ bacteriaMDS2uncsite Treatment +
Transpiration ~ Treatment #+ bacteriaMDS1uncsite #+ bacteriaMDS2uncsite

bacteriaMDS2uncsite~~0*FlowersperRosette+0*Transpiration
FlowersperRosette~~0*Transpiration 
LeafLengthmm~~0*Transpiration+0*bacteriaMDS2uncsite
'

#constrained by treatment, conditioned on site
mod.wmmm = '
bacteriaMDS1tresite ~ Treatment
#Biomassg ~ Treatment #+ bacteriaMDS1tresite
LeafLengthmm ~ Treatment# +bacteriaMDS1tresite
#LeafNumber ~ Treatment #+ bacteriaMDS1tresite
#FlowersperRosette ~ bacteriaMDS1tresite #Treatment +
Assimilation ~ Treatment +bacteriaMDS1tresite
Transpiration ~ Treatment + bacteriaMDS1tresite
'

mf.wmmm = sem(model = mod.wmmm, data = datEWMMMcstand) 
summary(mf.wmmm, fit.measures=TRUE,rsquare=T)
AIC(mf.wmmm)

#Direct: 
0.508+1.210+0.588
#Indirect: 
-0.736*0.966+-0.736*1.023

boxplot(datEWMMMcstand$Transpiration~as.factor(datEWMMMcstand$Treatment))
boxplot(datEWMMMb$LeafNumber~as.factor(datEWMMMb$Treatment))

ggplot(datEWMMMcstand, aes(x=bacteriaMDS1tresite, y=Transpiration,color=Treatment))+
  theme_classic()+
  geom_point()+
  geom_smooth(method="lm",se=F)

summary(lm(Transpiration~Treatment+bacteriaMDS1tresite+Site,data=datEWMMMcstand))
summary(lm(Transpiration~bacteriaMDS1tresite,data=datEWMMMcstand))










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

