#Data cleaning and merging


##### Lab data #####
labdat<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/Niwot_IndirectEffects_2021_LabData.csv",stringsAsFactors = T)
labdat2<-labdat%>%
  filter(SoilMoisturePercent>0)%>%
  filter(!(PlotID=="ES-FF2"&Replicate==1))%>%
  filter(!(PlotID=="AS-DM3"&Replicate==1))%>%
  filter(!(PlotID=="LE-DM1c"&Replicate==1))%>%
  filter(!(PlotID=="AE-MM1a"&Replicate==1))%>%
  filter(!(PlotID=="TC-WM1a"&Replicate==1))%>%
  filter(!(PlotID=="TE-WM3c"&Replicate==1))%>%
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
















##### Survey data #####

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

dat<-labdat2E%>%
  full_join(fielddatE)%>%
  full_join(licordat2E)%>%
  full_join(amf2E)
head(dat)
dim(dat)

#note eventually I want to reload these dataframes and have the site names listed out "Trough" and have Proj change to Treatment

