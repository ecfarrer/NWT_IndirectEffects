#Data cleaning and merging


##### Lab data #####
labdat<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/Niwot_IndirectEffects_2021_LabData.csv",stringsAsFactors = T)
labdat2<-labdat%>%
  filter(SoilMoisturePercent>0)%>%
  filter(!(PlotID=="ES_FF_2"&Replicate==1))%>%
  filter(!(PlotID=="AS_DM_3"&Replicate==1))%>%
  filter(!(PlotID=="LE_DM_1c"&Replicate==1))%>%
  filter(!(PlotID=="AE_MM_1a"&Replicate==1))%>%
  filter(!(PlotID=="TC_WM_1a"&Replicate==1))%>%
  filter(!(PlotID=="TE_WM_3c"&Replicate==1))%>%
  dplyr::select(PlotID:Replicate,SoilMoisturePercent,Biomassg)%>%
  group_by(PlotID,Site,PlotType,Treatment,SurveyAnalysis,ExperimentAnalysis,CommunityType,Chamber,Plot,PlantID)%>%
  summarise(across(SoilMoisturePercent:Biomassg,~mean(.x,na.rm=T)))%>%
  arrange(PlotID)
labdat2<-data.frame(labdat2)
labdat2$CommunityType<-factor(labdat2$CommunityType,levels=c("SB","WM","MM","DM","FF"))
View(labdat2)


##### Field data #####
fielddat<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/Niwot_IndirectEffects_2021_FieldData.csv",stringsAsFactors = T)
fielddat$CommunityType<-factor(fielddat$CommunityType,levels=c("SB","WM","MM","DM","FF"))

dim(labdat2)
dim(fielddat)



##### LICOR data #####

licordat<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/LiCordata.csv",stringsAsFactors = T)
head(licordat)

#The filtering of Acorrectedredo<50 removes two high FF numbers (TS_FF6, TS_FF2), one had a terrible picture (TS_FF6), the other had a lot of deschampsia in the frame (TS_FF2). i'm also removing AS_FF_1 b/c lots of des in frame.
licordat2<-licordat%>%
  dplyr::select(PlotID,CommunityType,LeafAreacm2:gtc)%>%  
  mutate(Ecorrected=Ecorrected*1000)%>%
  mutate(Ecorrectedredo=Ecorrectedredo*1000)%>%
  filter(!is.na(LeafAreacm2redo))%>% #if you want to look at stomatal conductance you need to do this, because there are weird negative numbers here
  filter(Acorrectedredo<50)%>% #Acorrected<40
  filter(PlotID!="AS_FF_1")

licordat2$CommunityType<-factor(licordat2$CommunityType,levels=c("SB","WM","MM","DM","FF"))

boxplot(licordat2$Acorrectedredo)

#duplicated(licordat2$SampleID)

#notes from annual report with Acorrect<40: TS-FF6 leaf area is small, AS-FF1 leaf area small, LS-FF1 is terrible pic, TS-FF2 has deschampsia in the chamber and is folded, EC-DM3b is folded
#filter(licordat2,Acorrected>40)

ggplot(licordat2,aes(y=Ecorrectedredo,x=CommunityType))+ #transpiration
  geom_boxplot()
ggplot(licordat2,aes(y=Acorrectedredo,x=CommunityType))+ #Assimilation
  geom_boxplot()
ggplot(licordat2,aes(y=gsw,x=CommunityType))+ #stomatal conductance to water
  geom_boxplot()



##### Root Staining #####

amf<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/amf2.csv",stringsAsFactors=T)
head(amf)
amf$CommunityType<-factor(amf$CommunityType,levels=c("SB","WM","MM","DM","FF"))

#If I'm really going to use this data, I need to go over some of these replicates, I just did a really quick determination of whether to average over reps or delete some. it was not at all clear what the write decision should be for a lot of the reps (i.e. there were no notes about why the rep was done and the results were different but not too different)
#I will write the reasoning just below if I keep all the reps or I'll write the reason after the filter() function if I filter some reps
filter(!(SampleID=="AC-MM3b"&Replicate==1))%>% #could delete first rep b/c arbuscules are a little low, but I'll keep it b/c it is not outside the range in the site/treatment
filter(!(SampleID=="ES-MM3"&Replicate==1))%>% #could delete first rep b/c vesicles are a little low but it is not extreme and it is within the range in the site/treatment
filter(!(SampleID=="LS-SB6"&Replicate==1))%>% #could delete first rep b/c everything is a bit low but not extreme and within the range in site/treatment
filter(!(SampleID=="TC-WM1c"&Replicate==1))%>% #the two reps are exactly the same
filter(!(SampleID=="TS-DM5"&Replicate==1))%>% #could delete b/c everything is a little low, but not outside the range, however the range is kind of oddly low in this site/treatment and many were done by Anna
filter(!(SampleID=="TS-SB1"&Replicate==1))%>% #could delete b/c arbuscules are a little low but not outside the range
filter(!(SampleID=="TS-SB2"&Replicate==1))%>% #could delete b/c arbuscules are a little low but not outside the range
  
amf2<-amf%>%
  filter(!(PlotID=="AC_MM_1c"&Replicate==2))%>% #comment that staining is weak and rep 2 has less colonization, so delete rep 2 b/c stain was probably even weaker than rep 1
  filter(!(PlotID=="AC_MM_3a"&Replicate==1))%>% #delete first two reps b/c arbuscules are really low
  filter(!(PlotID=="AC_MM_3a"&Replicate==2))%>% #delete first two reps b/c arbuscules are really low
  filter(!(PlotID=="AE_MM_1a"&Replicate==1))%>% #delete b/c kanley and total hyphae are really low
  filter(!(PlotID=="AE_MM_1b"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="AE_MM_1c"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="AE_MM_2a"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="AE_MM_2c"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="AS_MM_2"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="AS_MM_3"&Replicate==1))%>% #delete b/c kanley and total hyphae and arbuscules are really low
  filter(!(PlotID=="EE_DM_2b"&Replicate==1))%>% #delete b/c total hyphae are really low
  filter(!(PlotID=="EE_DM_3c"&Replicate==1))%>% #delete b/c total hyphae and arbuscules are really low
  filter(!(PlotID=="TE_WM_1a"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_1b"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_1c"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_2a"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_2b"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_2c"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_3a"&Replicate==1))%>% #delete b/c kanely and everything really low
  filter(!(PlotID=="TE_WM_3b"&Replicate==1))%>% #delete b/c kanely and arbuscules really low
  filter(!(PlotID=="TE_WM_3c"&Replicate==1))%>% #delete b/c kanely and everything really low
  #filter(Person!="KANELY")%>%
  #filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one. i took this out b/c there were so few WM. i'll take these out later anyway, not sure why this is here
  dplyr::select(PlotID,CommunityType:Replicate,Negative:HyphaeP)%>%
  group_by(PlotID,CommunityType,Plot,PlantID)%>%
  summarise(across(Negative:HyphaeP,~mean(.x,na.rm=T)))%>%
  arrange(PlotID)

dim(amf2)
#I'm missing one survey plot



##### Undergrad root staining #####

uamf<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/UndergradRootScoring.csv",stringsAsFactors=T)
head(uamf)

#This has not been cleaned, there are some notes about weak stains that I could address and delete. I just want to look at the first pass patterns

#amf$CommunityType<-factor(amf$CommunityType,levels=c("SB","WM","MM","DM","FF"))
uamf2<-uamf%>%
  filter(PlotID!="TC_MM_2")%>% #this might really be TS_MM_2 but for now I'll just delete it
  dplyr::select(PlotID,Rep,Negative:HyphaeP)%>%
  group_by(PlotID)%>%
  summarise(across(Negative:HyphaeP,~mean(.x,na.rm=T)))%>%
  arrange(PlotID)%>%
  left_join(labdat2)
uamf2$SiteTreatment<-paste(uamf2$Site,uamf2$Treatment,sep="_")
uamf2$SiteTreatment<-factor(uamf2$SiteTreatment,levels=c("Trough_Control","Trough_Experimental","Trough_","Audubon_Control","Audubon_Experimental","Audubon_","Lefty_Control","Lefty_Experimental","Lefty_","EastKnoll_Control","EastKnoll_Experimental","EastKnoll_"))
uamf2$Site<-factor(uamf2$Site,levels=c("Trough","Audubon","Lefty","EastKnoll"))


##### Merge all datasets #####

dat<-labdat2%>%
  full_join(fielddat)%>%
  full_join(licordat2)%>%
  full_join(amf2)
dat$Site<-factor(dat$Site,levels=c("Trough","Audubon","Lefty","EastKnoll"))
dat$SiteCommunityType<-paste(dat$Site,dat$CommunityType,sep="")
dat$SiteCommunityType<-factor(dat$SiteCommunityType,levels=c("TroughSB","TroughWM","TroughMM","TroughDM","TroughFF","AudubonSB","AudubonMM","AudubonDM","AudubonFF","LeftySB","LeftyMM","LeftyDM","LeftyFF","EastKnollSB","EastKnollMM","EastKnollDM","EastKnollFF"))
head(dat)
dim(dat)

#I caught and fixed one copy and paste error from the field data leaf area into the licor data leaf area. and I checked that all the NAs/differences between the field data leaf areas and the licor data leaf areas are accounted for. the ones that are NAs in the licor data leaf area are samples that were removed b/c they had suspicious values for Acorrectedredo or had a deschampsia filling most of the frame so suspect
dat$LeafAreacm2redo-dat$RedoLeafAreacm2fielddata
dat[143:147,c("PlotID","LeafAreacm2redo")]
dat[143:147,c("PlotID","RedoLeafAreacm2fielddata")]



##### Adding some columns to microbe data sets #####
sample_data(datITSS5c)$SiteCommunityType<-paste(sample_data(datITSS5c)$Site,sample_data(datITSS5c)$CommunityType,sep="")
sample_data(datITSS5c)$SiteCommunityType<-factor(sample_data(datITSS5c)$SiteCommunityType,levels=c("AudubonSB","AudubonMM","AudubonDM","AudubonFF","EastKnollSB","EastKnollMM","EastKnollDM","EastKnollFF","LeftySB","LeftyMM","LeftyDM","LeftyFF","TroughSB","TroughWM","TroughMM","TroughDM","TroughFF"))

sample_data(datITSS5c)$SiteTreatment<-paste(sample_data(datITSS5c)$Site,sample_data(datITSS5c)$Treatment,sep="")
sample_data(datITSS5c)$SiteTreatment<-factor(sample_data(datITSS5c)$SiteTreatment,levels=c("AudubonControl","AudubonExperimental","EastKnollControl","EastKnollExperimental","LeftyControl","LeftyExperimental","TroughControl","TroughExperimental","Audubon","EastKnoll","Lefty","Trough"))

sample_data(dat16SS5c)$SiteCommunityType<-paste(sample_data(dat16SS5c)$Site,sample_data(dat16SS5c)$CommunityType,sep="")
sample_data(dat16SS5c)$SiteCommunityType<-factor(sample_data(dat16SS5c)$SiteCommunityType,levels=c("AudubonSB","AudubonMM","AudubonDM","AudubonFF","EastKnollSB","EastKnollMM","EastKnollDM","EastKnollFF","LeftySB","LeftyMM","LeftyDM","LeftyFF","TroughSB","TroughWM","TroughMM","TroughDM","TroughFF"))

sample_data(dat16SS5c)$SiteTreatment<-paste(sample_data(dat16SS5c)$Site,sample_data(dat16SS5c)$Treatment,sep="")
sample_data(dat16SS5c)$SiteTreatment<-factor(sample_data(dat16SS5c)$SiteTreatment,levels=c("AudubonControl","AudubonExperimental","EastKnollControl","EastKnollExperimental","LeftyControl","LeftyExperimental","TroughControl","TroughExperimental","Audubon","EastKnoll","Lefty","Trough"))

sample_data(dat16SS5c)$CommunityType<-factor(sample_data(dat16SS5c)$CommunityType,levels=c("SB","WM","MM","DM","FF"))




##### Getting Survey data together #####

#think about how I want to do the subsetting for the survey. do I want to include all the experimental controls or just one per plot. As of now I will include only the a's from the control plots for all sites except the WM in Trough fow which I will use all of them

#In general we only have WM from the control (and experimental plots), not from any plots outside the turf area
#For lab/field data, we are missing licor from TC_WM_1b, missing amf from ES_DM_3
#For fungi we are missing TC_WM_1a roots,TC_WM_1b roots, TS_FF_5 soil
#For bacteria we are missing AS_SB_1 roots, LS_FF_3 roots, LS_SB_5 roots, TC_WM_2c roots

datS<-dat%>%
  filter(SurveyAnalysis=="Survey")%>%
  #filter(PlantID%nin%c("1b","1c","2b","2c","3b","3c"))%>% #take out samples from the same control plot, only use one
  arrange(Site,CommunityType,Treatment,PlantID)
#105 samples

# datS%>%
#   filter(Treatment=="Control")
# 
# datS%>%
#   group_by(Site,CommunityType)%>%
#   tally()
# 
# temp<-data.frame(sample_data(dat16SS5c))%>%
#   filter(SurveyAnalysis=="Survey")%>%
#   arrange(Site,CommunityType,SampleType)%>%
#   #filter(!(Rep%in%c("1b","1c","2b","2c","3b","3c")))%>%
#   group_by(Site,CommunityType,SampleType)%>%
#   tally()
# as.data.frame(temp)


datITSS5cS<-datITSS5c%>%
  subset_samples(SurveyAnalysis=="Survey")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
#207 samples, missing 3

dat16SS5cS<-dat16SS5c%>%
  subset_samples(SurveyAnalysis=="Survey")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
#206 samples, missing 4





##### Getting Experiment data together #####

#Lab/field, missing licor from TC_WM_1b
#Fungi,
#Bacteria,


datE<-dat%>%
  filter(ExperimentAnalysis=="Experiment")%>%
  arrange(Site,CommunityType,Treatment,PlantID)%>%
  mutate(CommunityTreatment=factor(paste(CommunityType,Treatment,sep="_"),levels=c("WM_Control","WM_Experimental","MM_Control","MM_Experimental","DM_Control","DM_Experimental")))%>%
  mutate(SiteTreatment=factor(paste(Site,Treatment,sep="_"),levels=c("Trough_Control","Trough_Experimental","Audubon_Control","Audubon_Experimental","Lefty_Control","Lefty_Experimental","EastKnoll_Control","EastKnoll_Experimental")))%>%
  mutate(SitePlot=factor(paste(Site,Plot,sep="_")))
datE$MoistureType=datE$CommunityType
ind<-datE$MoistureType=="WM"
datE$MoistureType[ind]<-"MM"
datE$MoistureTreatment<-paste(datE$MoistureType,datE$Treatment,sep="_")
dim(datE)
#72 sample (=9*2*4)


datITSS5cE<-datITSS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
#140 samples, missing 4, EE_DM_1c soil, TC_WM_1a roots,TC_WM_1b roots, TE_WM_3c roots


dat16SS5cE<-dat16SS5c%>%
  subset_samples(ExperimentAnalysis=="Experiment")%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
#140 samples, missing 4, AC_MM_1c roots, AE_MM_3b soil, AC_MM_1b soil, TC_WM_2c roots


# temp<-data.frame(sample_data(dat16SS5cE))%>%
#   group_by(Site,CommunityType,SampleType)%>%
#   tally()
# as.data.frame(temp)




##### Get trough survey and experiment data together #####

datITSS5cTrough<-datITSS5c%>%
  subset_samples(Site=="Trough")%>%
  subset_samples(CommunityType%in%c("WM","MM"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sample_data(datITSS5cTrough)$Group<-sample_data(datITSS5cTrough)$Treatment
ind<-sample_data(datITSS5cTrough)$Group==""
sample_data(datITSS5cTrough)$Group[ind]<-"Survey"

dat16SS5cTrough<-dat16SS5c%>%
  subset_samples(Site=="Trough")%>%
  subset_samples(CommunityType%in%c("WM","MM"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sample_data(dat16SS5cTrough)$Group<-sample_data(dat16SS5cTrough)$Treatment
ind<-sample_data(dat16SS5cTrough)$Group==""
sample_data(dat16SS5cTrough)$Group[ind]<-"Survey"

#Audubon, experimental are all moist meadow, also have 3 more MM in survey i should delete
datITSS5cAudubon<-datITSS5c%>%
  subset_samples(Site=="Audubon")%>%
  subset_samples(CommunityType%in%c("MM","DM"))%>%
  subset_samples(!(CommunityType%in%c("MM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sample_data(datITSS5cAudubon)$Group<-sample_data(datITSS5cAudubon)$Treatment
ind<-sample_data(datITSS5cAudubon)$Group==""
sample_data(datITSS5cAudubon)$Group[ind]<-"SurveyDM"

dat16SS5cAudubon<-dat16SS5c%>%
  subset_samples(Site=="Audubon")%>%
  subset_samples(CommunityType%in%c("MM","DM"))%>%
  subset_samples(!(CommunityType%in%c("MM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sample_data(dat16SS5cAudubon)$Group<-sample_data(dat16SS5cAudubon)$Treatment
ind<-sample_data(dat16SS5cAudubon)$Group==""
sample_data(dat16SS5cAudubon)$Group[ind]<-"SurveyDM"

#Lefty, experimental are all DM, also have 3 more DM in survey i should delete
datITSS5cLefty<-datITSS5c%>%
  subset_samples(Site=="Lefty")%>%
  subset_samples(CommunityType%in%c("DM","FF"))%>%
  subset_samples(!(CommunityType%in%c("DM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
#sample_data(datITSS5cAudubon)$Group<-sample_data(datITSS5cAudubon)$Treatment
#ind<-sample_data(datITSS5cAudubon)$Group==""
#sample_data(datITSS5cAudubon)$Group[ind]<-"SurveyDM"

dat16SS5cLefty<-dat16SS5c%>%
  subset_samples(Site=="Lefty")%>%
  subset_samples(CommunityType%in%c("DM","FF"))%>%
  subset_samples(!(CommunityType%in%c("DM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
# sample_data(dat16SS5cAudubon)$Group<-sample_data(dat16SS5cAudubon)$Treatment
# ind<-sample_data(dat16SS5cAudubon)$Group==""
# sample_data(dat16SS5cAudubon)$Group[ind]<-"SurveyDM"

#EastKnoll, experimental are all DM, also have 3 more DM in survey i should delete
datITSS5cEastKnoll<-datITSS5c%>%
  subset_samples(Site=="EastKnoll")%>%
  subset_samples(CommunityType%in%c("DM","FF"))%>%
  subset_samples(!(CommunityType%in%c("DM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

dat16SS5cEastKnoll<-dat16SS5c%>%
  subset_samples(Site=="EastKnoll")%>%
  subset_samples(CommunityType%in%c("DM","FF"))%>%
  subset_samples(!(CommunityType%in%c("DM")&PlotType=="Survey"))%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

