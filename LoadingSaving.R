
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
library(NetCoMi) #networks
#library(deming)#for deming regression with x and y error
#library(emmeans)#for least squares means
#library(ggforce)#for some ellipses on ordinations
#library(lavaan)

options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")


#This is my current analysis workspace
save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis.RData")
save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis2.RData")
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspaceanalysis2.RData")

#This is for my first runthrough only with roots for the annual report
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspace1analysis.RData")

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


##### Network analysis loading and saving nonsense #####

#I tried this with R version 4.3.3 - but it didn't work. NetCoMi still needs metagenomeSeq version 1.47 and doing it with R version 4.3.3 and bioconductor 3.18 installs metagenomeSeq version 1.43. So I think I need to use R 4.4 and just wait until the metagenomeSeq installer works properly

#metagenomeSeq@1.47.0 should be 1.47.0, i currently have 1.46
#there is something wrong with the biocmanager installer for the current version of metagenomeSeq

# Now I'm trying with the updated R 4.4.2, and installing biocmanager 3.20. still the metagenomeSeq that can be installed is 1.46.0  This is messed up.

sessionInfo()
packageVersion("metagenomeSeq") 

#gfortran needs to be downloaded from the r website: https://cran.r-project.org/bin/macosx/tools/
#Bioconductor version needs to match R version. 3.18 matches with R 4.3.3

install.packages("devtools")
install.packages("BiocManager")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")

#install.packages("BiocManager", repos = "https://cloud.r-project.org")

#install.packages("devtools")

devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")

devtools::install_github("stefpeschel/NetCoMi", 
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

installNetCoMiPacks()

# BiocManager::install("metagenomeSeq")
# BiocManager::install("metagenomeSeq",force=TRUE)
# packageVersion("metagenomeSeq") #this installed 1.46.0

#BiocManager::install(version = "3.19", force=T)
#BiocManager::install("metagenomeSeq",force = T)#,version="3.19"

#devtools::install_github("Bioconductor-mirror/metagenomeSeq") #doesn't work
# devtools::install_github("Bioc/metagenomeSeq") #installs version 1.49
# devtools::install_github("Bioc/metagenomeSeq@RELEASE_3_20") #installs version 1.48
# devtools::install_github("Bioc/metagenomeSeq@RELEASE_3_19") #installs version 1.46
# devtools::install_github("HCBravoLab/metagenomeSeq") #installs version 1.31

# devtools::install_github("stefpeschel/NetCoMi", 
#                          dependencies = c("Depends", "Imports", "LinkingTo"),
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))
library(NetCoMi)

# library(SpiecEasi)
# library(SPRING)
# library(metagenomeSeq)
# library(mixedCCA)




