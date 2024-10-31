
library(phyloseq)
library(dada2)
library(tidyverse)
library(vegan)
library(permute)
library(plotrix)
#library(fossil) #for sorenson()
library(nlme)
library(cowplot)
#library(car)
library(MASS)
#library(pscl)
library(multcomp)
library(Hmisc)
#library(deming)#for deming regression with x and y error
library(emmeans)#for least squares means
#library(ggforce)#for some ellipses on ordinations
library(lavaan)

options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

#This is for my first runthrough only with roots for the annual report
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspace1analysis.RData")

#This is my current analysis workspace
save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis.RData")
save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis2.RData")
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis2.RData")


#This is updated for fungal and bacterial bioinformatics
save.image(file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspacecleandata.Rdata")
 


#Bioinformatics on fungi 
workspacebioinformaticsITS.Rdata
#Bioinformatics on bacteria
workspacebioinformatics16S.RData

#After bioinformatics only necessary dataframes
workspacecleandata.Rdata

rm(list=setdiff(ls(), c("datITSS3","datITSS4","datITSS5", "datITSS5c","datITSS5otu","datITSS5cotu","funguild2")))
datITSS3 #all samples just after decontam
datITSS4 #<4000 reads removed
datITSS5 # rarefied rel abun
datITSS5c # rarefied counts
datITSS5otu #dataframe otu table rel abun
datITSS5cotu  #dataframe otu table counts
funguild2

rm(list=setdiff(ls(),c("dat16SS3a","dat16SS4","dat16SS5","dat16SS5c","dat16SS4_r4903","dat16SS5_r4903","dat16SS5c_r4903","dat16SS5otu","dat16SS5cotu","dat16SS5otu_r4903","dat16SS5cotu_r4903","faprotaxOTU","faprotaxSample","faprotaxSamplelong")))
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspacecleandata.Rdata")
save.image(file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspacecleandata.Rdata")
dat16SS3a #all samples after decontam and filtering for just archaea and bacteria
dat16SS4 #<2000 reads removed
dat16SS5 # rarefied rel abun
dat16SS5c # rarefied counts
dat16SS4_r4903 #<4900 reads removed
dat16SS5_r4903 # rarefied rel abun
dat16SS5c_r4903 # rarefied counts
dat16SS5otu
dat16SS5cotu
dat16SS5otu_r4903
dat16SS5cotu_r4903
faprotaxOTU #The raw output of what functions the different OTUs have
faprotaxSample#The abundance weighted percent of each function in each sample
faprotaxSamplelong#The abundance weighted percent of each function in each sample in long format with categories (n cycling, c cycling)


#Analysis workspaces 
workspace5.Rdata

#OLD BELOW


