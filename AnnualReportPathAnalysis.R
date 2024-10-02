#Path analysis


##### Experiment WM models #####

datEWM<-dat%>%
  filter(Proj!="S",Community=="WM")
temp<-sample_data(tempphyEWM)
temp$SampleID<-gsub("_", "-", temp$PlotID)
temp$SampleID<-sub("(-.*?)-", "\\1", temp$SampleID)

datEWMb<-datEWM%>%
  full_join(temp,by="SampleID")%>%
  filter(!is.na(MDS1))%>%
  filter(!is.na(Acorrected))
datEWMb$Proj<-as.numeric(as.character(recode_factor(datEWMb$Proj,"C"=0,"E"=1)))

#do i need to pre-standardize everything? yes it looks like particularly due to the Proj binary variable which is not standardized or centered (?) using the std.ov=T
datEWMc<-datEWMb%>%
  dplyr::select(Proj, NMDS1, NMDS2, MDS1, MDS2, Biomassg,LeafNumber,LeafLengthmm, FlowersperRosette, Ecorrected, Acorrected)
  
datEWMcstand<-apply(datEWMc,2,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})


#Using NMDS, I like NMDS b/c the points are not as far apart (like outliers) on the ordination axis (which makes some outlier points have too much leverage). I took out A corrected at the end b/c I realized that only growth/flowers was originally going to be in the path diagram
mod.wm = '
NMDS1 ~ Proj
#NMDS2 ~ Proj
Biomassg ~ Proj + NMDS1 #+ NMDS2  
LeafLengthmm ~ Proj + NMDS1 #+ NMDS2 
LeafNumber ~  NMDS2#+ NMDS1 Proj  +
FlowersperRosette ~ Proj #+ NMDS1# + NMDS2
#Acorrected ~ Proj + NMDS1 + NMDS2
#Ecorrected ~ Proj # + NMDS2#+ NMDS1
'

mod.est.wm = sem(model = mod.wm, std.ov=T, data = datEWMcstand) #datEWMb
AIC(mod.est.wm)
summary(mod.est.wm, fit.measures=TRUE,rsquare=T)
AIC(mod.est.wm)

#Direct: 
0.508+1.210+0.588
#Indirect: 
-0.736*0.966+-0.736*1.023





#Using NMDS, removed some covariances. doesn't change much at all
mod.wm = '
NMDS1 ~ Proj
#NMDS2 ~ Proj
Biomassg ~ Proj + NMDS1 # + NMDS2
LeafLengthmm ~ Proj + NMDS1 # + NMDS2 
LeafNumber ~  + NMDS2 #+ NMDS1 Proj 
FlowersperRosette ~ Proj # + NMDS1 # + NMDS2
Acorrected ~ Proj + NMDS1 + NMDS2
#Ecorrected ~ Proj # + NMDS2 #+ NMDS1
Biomassg~~0*Acorrected#+0*Ecorrected
LeafLengthmm~~0*Acorrected#+0*Ecorrected
LeafNumber~~0*Acorrected#+0*Ecorrected
FlowersperRosette~~0*Acorrected#+0*Ecorrected
'
mod.est.wm = sem(model = mod.wm, std.ov=T, data = datEWMb)
AIC(mod.est.wm)
summary(mod.est.wm, fit.measures=TRUE,rsquare=T)
AIC(mod.est.wm)

#Using MDS
mod.wm = '
MDS1 ~ Proj
MDS2 ~ Proj
#Biomassg ~ Proj  #+ MDS2#+ MDS1
#LeafLengthmm ~ Proj #+ MDS1 #+ MDS2 
#LeafNumber ~ MDS2 #+ MDS1# Proj  +
FlowersperRosette ~ Proj + MDS1 #+ MDS2
Acorrected ~MDS2 # Proj + MDS1 + 
Ecorrected ~   MDS2#Proj +MDS1 +
'

mod.est.wm = sem(model = mod.wm, std.ov=T, data = datEWMb)
AIC(mod.est.wm)
summary(mod.est.wm, fit.measures=TRUE,rsquare=T)
AIC(mod.est.wm)



plot(datEWMb$MDS1,datEWMb$FlowersperRosette)
boxplot(datEWMb$FlowersperRosette~datEWMb$Proj)

ggplot(datEWMb,aes(x=NMDS1,y=Biomassg))+
  geom_point()+
  geom_smooth(method="lm")



# library(semPlot)
# semPaths(
#   object = mod.est.wm,
#   what = "path",
#   whatLabels = "par",style="ram",
#   rotation=2
# )



mod.wm = '
MDS1 ~ Proj
MDS2 ~ Proj
Biomassg ~ Proj + MDS1 + MDS2
LeafLengthmm ~ Proj + MDS1 + MDS2 
LeafNumber ~ Proj + MDS1 + MDS2
FlowersperRosette ~ Proj + MDS1 + MDS2
Acorrected ~ Proj + MDS1 + MDS2
Ecorrected ~ Proj + MDS1 + MDS2
'

mod.wm = '
NMDS1 ~ Proj
NMDS2 ~ Proj
Biomassg ~ Proj + NMDS1 + NMDS2
LeafLengthmm ~ Proj + NMDS1 + NMDS2 
LeafNumber ~ Proj + NMDS1 + NMDS2
FlowersperRosette ~ Proj + NMDS1 + NMDS2
Acorrected ~ Proj + NMDS1 + NMDS2
Ecorrected ~ Proj + NMDS1 + NMDS2
'







##### Experiment DM models #####

datEDM<-dat%>%
  filter(Proj!="S",Community=="DM")
temp<-sample_data(tempphyEDM)
temp$SampleID<-gsub("_", "-", temp$PlotID)
temp$SampleID<-sub("(-.*?)-", "\\1", temp$SampleID)

datEDMb<-datEDM%>%
  full_join(temp,by="SampleID")%>%
  filter(!is.na(MDS1))%>%
  filter(!is.na(Acorrected))
datEDMb$Proj<-as.numeric(as.character(recode_factor(datEDMb$Proj,"C"=0,"E"=1)))

#do i need to pre-standardize everything? yes it looks like particularly due to the Proj binary variable which is not standardized or centered (?) using the std.ov=T
datEDMc<-datEDMb%>%
  dplyr::select(Proj, NMDS1, NMDS2, MDS1, MDS2, Biomassg,LeafNumber,LeafLengthmm, FlowersperRosette, Ecorrected, Acorrected)

datEDMcstand<-apply(datEDMc,2,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})


#Using NMDS, I like NMDS b/c the points are not as far apart (like outliers) on the ordination axis (which makes some outlier points have too much leverage). I took out A corrected at the end b/c I realized that only growth/flowers was originally going to be in the path diagram
mod.dm = '
#NMDS1 ~ Proj
NMDS2 ~ Proj
Biomassg ~ Proj +NMDS1# + NMDS2   
#LeafLengthmm ~ Proj  #+ NMDS1 + NMDS2
LeafNumber ~   NMDS1 #+ NMDS2  Proj +
FlowersperRosette ~ Proj #+ NMDS2 #+ NMDS1 
Acorrected ~ Proj  + NMDS2#+ NMDS1
Ecorrected ~  NMDS1 #+ NMDS2 #Proj +
'

mod.est.dm = sem(model = mod.dm, std.ov=T, data = datEDMcstand) 
AIC(mod.est.dm)
summary(mod.est.dm, fit.measures=TRUE,rsquare=T)
AIC(mod.est.dm)

#Direct: 
0.387+0.325
#Indirect: 
0.446*-0.516


boxplot(datEDMb$FlowersperRosette~datEDMb$Proj)

ggplot(datEDMc,aes(x=NMDS1,y=Biomassg))+
  geom_point()+
  geom_smooth(method="lm")



