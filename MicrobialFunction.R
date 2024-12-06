#Function

##### Fungal Traits #####
funtraits2
colnames(funtraits2)

m1<-funtraits2%>%
  filter(PlotType!="Survey")%>%
  group_by(SampleType,CommunityType,Treatment)%>%
  summarise(mean=mean(PrimandSecplant_pathogenFT), se=std.error(PrimandSecplant_pathogenFT))
m1
m1$CommunityType<-factor(m1$CommunityType,levels=c("WM","MM","DM"))
#m1$Site<-factor(m1$Site,levels=c("Trough","Audubon","Lefty","EastKnoll"))
#m1$CommunityType<-recode_factor(m1$CommunityType,"C"="Control","E"="Treatment")

ggplot(data=m1, aes(x=Treatment, y=mean,color=SampleType))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("")+
  #  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "right",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(CommunityType),strip.position = "bottom",nrow=1)

#do a loop for effect of moisture * treatment on each dependent variable with site/row as random effect


#for roots columns 16:142
funtraits2E<-funtraits2%>%
  filter(PlotType!="Survey")%>%
  filter(SampleType=="roots")
ind<-which(colSums(funtraits2E[,16:156])==0)
funtraits2E[,ind+15]<-NULL
#only necessary for Treatment*Community type interaction
funtraits2E$Animal_biotrophic_capacity.arthropodassociated<-NULL
colnames(funtraits2E)

output<-data.frame(var=rep(NA,126),treat=rep(NA,126),moisttype=rep(NA,126),interaction=rep(NA,126))
for(i in 16:141){
  tempvar<-colnames(funtraits2E)[i]
  tempdat<-funtraits2E[,i]
  m0<-lme(tempdat ~ Treatment*CommunityType, random=~1|Site/Plot/Chamber, na.action=na.omit,data = funtraits2E)
  output[i-15,1]<-tempvar
  output[i-15,2:4]<-anova(m0,type="margin")[2:4,"p-value"]
}
View(output)

#for soil columns 16:145
funtraits2E<-funtraits2%>%
  filter(PlotType!="Survey")%>%
  filter(SampleType=="soil")
ind<-which(colSums(funtraits2E[,16:156])==0)
funtraits2E[,ind+15]<-NULL
#for moisturetype, not community type
#funtraits2E$Ectomycorrhiza_exploration_type.shortdistance_coarse<-NULL
colnames(funtraits2E)

output<-data.frame(var=rep(NA,131),treat=rep(NA,131),moisttype=rep(NA,131),interaction=rep(NA,131))
for(i in 16:146){
  tempvar<-colnames(funtraits2E)[i]
  tempdat<-funtraits2E[,i]
  m0<-lme(tempdat ~ Treatment*CommunityType, random=~1|Site/Plot/Chamber, na.action=na.omit,data = funtraits2E)
  output[i-15,1]<-tempvar
  output[i-15,2:4]<-anova(m0,type="margin")[2:4,"p-value"]
}
View(output)
cbind(m1$SampleNumber,m1$PlotID,m1$Plot,m1$Chamber)





##### Bacteria FAPROTAX #####
bactraits2
colnames(bactraits2)

m1<-bactraits2%>%
  filter(PlotType!="Survey")%>%
  group_by(SampleType,Site,Treatment)%>%
  summarise(mean=mean(anaerobic_chemoheterotrophy), se=std.error(anaerobic_chemoheterotrophy))
#m1$CommunityType<-factor(m1$CommunityType,levels=c("WM","MM","DM"))
m1$Site<-factor(m1$Site,levels=c("Trough","Audubon","Lefty","EastKnoll"))

#Community type
ggplot(data=m1, aes(x=Treatment, y=mean,color=SampleType))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("")+
  #  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "right",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(CommunityType),strip.position = "bottom",nrow=1)

#Site
ggplot(data=m1, aes(x=Treatment, y=mean,color=SampleType))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("")+
  #  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "right",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(Site),strip.position = "bottom",nrow=1)

#do a loop for effect of moisture * treatment on each dependent variable with site/row as random effect


#for roots columns 16:80
colnames(bactraits2E)
bactraits2E<-bactraits2%>%
  filter(PlotType!="Survey")%>%
  filter(SampleType=="roots")
ind<-which(colSums(bactraits2E[,16:80])==0)
bactraits2E[,ind+15]<-NULL
#only necessary for Treatment*Community type interaction
#funtraits2E$Animal_biotrophic_capacity.arthropodassociated<-NULL
colnames(bactraits2E)

output<-data.frame(var=rep(NA,56),treat=rep(NA,56),moisttype=rep(NA,56),interaction=rep(NA,56))
for(i in 16:71){
  tempvar<-colnames(bactraits2E)[i]
  tempdat<-bactraits2E[,i]
  m0<-lme(tempdat ~ Treatment*CommunityType, random=~1|Site/Plot/Chamber, na.action=na.omit,data = bactraits2E)
  output[i-15,1]<-tempvar
  output[i-15,2:4]<-anova(m0,type="margin")[2:4,"p-value"]
}
View(output)


#for soil columns 16:145
bactraits2E<-bactraits2%>%
  filter(PlotType!="Survey")%>%
  filter(SampleType=="soil")
ind<-which(colSums(bactraits2E[,16:80])==0)
bactraits2E[,ind+15]<-NULL
#for moisturetype, not community type
#funtraits2E$Ectomycorrhiza_exploration_type.shortdistance_coarse<-NULL
colnames(bactraits2E)

output<-data.frame(var=rep(NA,62),treat=rep(NA,62),moisttype=rep(NA,62),interaction=rep(NA,62))
for(i in 16:77){
  tempvar<-colnames(bactraits2E)[i]
  tempdat<-bactraits2E[,i]
  #m0<-lme(tempdat ~ Treatment*CommunityType, random=~1|Site/Plot/Chamber, na.action=na.omit, data = bactraits2E)
  m0<-lme(tempdat ~ Treatment*Site, random=~1|SitePlot/Chamber, na.action=na.omit, data = bactraits2E)
  output[i-15,1]<-tempvar
  output[i-15,2:4]<-anova(m0,type="margin")[2:4,"p-value"]
}
View(output)

sort(output$moisttype)
sort(p.adjust(output$treat, method="fdr"))

