#2021 bioiinformatics for Root and soil bacteria 16S samples for NSF Eager Indirect Effects Project

#Be sure to go to terminal and run the following to open the R version 4.4.1
rig default 4.4-x86_64

#First filter and infer sequence variants for each set: Roots, Soil
#For each set, run through fist part of script, CHANGE DIRECTORY folder

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspacebioinformatics16S.RData")
save.image(file = "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/workspacebioinformatics16S.RData")

# BiocManager::install(version = '3.16') #old
# BiocManager::version() #the one for R 4.4 is 3.19

#BiocManager::install(("decontam"))
#BiocManager::install(("phyloseq"))

#newest dada2 found here: https://bioconductor.org/packages/release/bioc/html/dada2.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

packageVersion("dada2") 
# Should be 1.32.0 #old was 1.26.0

library(dada2)
library(ShortRead)
library(Biostrings)
library("phyloseq")
library("plyr")
library("ggplot2")
library(tidyverse)
library(nlme)
library(vegan)
#library(reshape)
#library(stringr)
library(plotrix)
#library(data.table)
library(decontam)

#a pipeline to follow for 16S and dada2 in R
https://shibalytics.com/teaching/16s/microbiome_analysis_part_i
https://shibalytics.com/teaching/16s/microbiome_analysis_part_ii



#in unix add r or s to beginning of file names
for i in *; do mv "$i" r"$i"; done;ls -l
for i in *; do mv "$i" s"$i"; done;ls -l

#change path for each sequence set
path<- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/Roots16S"
path<- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/Soil16S" 

list.files(path)

#Look at a sequence from one sample
# test <- readFastq("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/Roots16S/r1_S1_L001_R1_001.fastq.gz")
# test
# width(test)
# str(sread(test)[1])
# as.character(sread(test)[1])


#separate forward and reverse
fnFs <- sort(list.files(path,pattern= "_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path,pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

#look at quality here?? I feel like it is harder to see what is happeneing compared to after you trim primers off
plotQualityProfile(fnFs[1:6])+ 
  geom_hline(yintercept=30)+
  geom_vline(xintercept=280)+#261+19
  geom_vline(xintercept=273)+#254+19
  geom_vline(xintercept=244)#300-56
plotQualityProfile(fnRs[1:6])+
  geom_hline(yintercept=30)+
  geom_vline(xintercept=225)#300-75

#Identify Primers
#forward primer: GTGYCAGCMGCCGCGGTAA
#reverse primer: GGACTACNVGGGTWTCTAAT

#515F + adapt
#CACTCTTTCCCTACACGACGCTCTTTCGATCTGTGYCAGCMGCCGCGGTAA
#806R + adapt
#GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGACTACNVGGGTWTCTAAT

#amplicon should be 390bp based on earth microbiome prooject. but based on our previous runs the final amplicons are 253 (maybe b/c they don't have primers/adapters, but I didn't check)

FWD <- "GTGYCAGCMGCCGCGGTAA" 
REV <- "GGACTACNVGGGTWTCTAAT" 

#verify that we have the right primers and correct orientation 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


#The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to "pre-filter" the sequences just to remove those with Ns, but perform no other filtering. I trimmed the ends of forward and reverse reads so that the median quality score was roughly 30.
#start 7:12pm, end 7:19pm
fnFs.filtN <- file.path(path,"filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path,"filtN", basename(fnRs))

#Initial, no trimming
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#For roots
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, trimRight = c(28,75))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, trimRight = c(56,75)) #try more aggressive trimming due b/c we are getting 254 rather than 253 length, that is odd. It is becaue there are 301 bp not 300, but in any case If I'm trimming more than abou 10 bp off the end, I should over trim so that I'm not trimming only part of the primer. Because if I leave 1-2 bp of primer after trimming, it won't get recongized as primer and it won't get trimmed. I will use this!

#For soil
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, trimRight = c(51,80))


#Count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient, assuming all the files were created using the same library preparation, so we'll just process the first sample.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Remove primers
cutadapt <- "/Users/farrer/opt/miniconda3/envs/qiime2-2023.2/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R

#Create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command
#critical parameters are the primers, and they need to be in the right orientation,
# i.e. the FWD primer should have been matching the forward-reads in its forward orientation, and the REV primer should have been matching the reverse-reads in its forward orientation. 
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt, on 4 cores 4:15-4:20
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 20, "-j", 4,# -n 2 required to remove FWD and REV from reads, -m means it deletes reads of less than 20 length, -j is number of cores (i didn't have this flag before, trying it out to see if it is faster)
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#sanity check:count presence of primers in cut-adapted files
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#green is mean, median is solid peach (go with median >=30)
#note: I just realized that the primer is taken off, so the forward primer is 19 and cutadapt took off 19 from all the forward reads so the read length is 281 after the untrimmed go through. so I should get the x value and then do 281-x, not 300-x. and the reverse primer is 20bp, so 280-x (really should be 301-19=282 and 301-20=281)
plotQualityProfile(cutFs[31:36])+ 
  geom_vline(xintercept=281)+ #281- taking into consideration the forward primer
  geom_hline(yintercept=30)+
  geom_vline(xintercept=261)+ #taking into consideration the read through reverse primer
  geom_vline(xintercept=254)+ #the amplicon should be mostly 253
  #geom_vline(xintercept=225)#trying more aggressive trimming roots
  geom_vline(xintercept=230) #soils
#for roots try trim at 281-253=28, then try trimming at 281-225=56
#for soil, trim at 282-230=52, do 51 just to make it the same total amount trimmed as the roots

plotQualityProfile(cutRs[1:6])+
  geom_vline(xintercept=280)+ #280-
  geom_hline(yintercept=30)+
  geom_vline(xintercept=261)+
  geom_vline(xintercept=254)+
  #geom_vline(xintercept=205)#roots
  geom_vline(xintercept=200)#soils

#for roots trim, due to low quality, at 280-205=75
#for soil, trim at 281-200=81, do 80 to make it the same total amount trimmed as the roots

#filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


#Standard filtering parameters maxN=0 (already filtered), truncQ=2, rm.phix=TRUE, maxEE=TRUE. 
#enforce min length of 50 bp
#start 10pm, end 10:06
#trimming low quality bases off at the first step results in over 2* the number of reads!
#I get about 1000 more reads per samples when using 85,75, but below after joining paired reads it is more variable
outroots2875 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
outroots5675 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
outsoil5180 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(outroots2875)
head(outroots5675)
head(outsoil5180)
cbind(outroots4447[,2],outroots8575[,2])
#write.csv(outroots, "readsout.csv")

#inspect read quality profiles
plotQualityProfile(filtFs[70])
plotQualityProfile(filtRs[70])

#learn Error rates, start 10:20, end 10:22
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference, 4 min
#At this step the core sample inference algorithm is applied to the dereplicated data
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge paired reads, 2 min
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) #default minOverlap=12

#construct sequence table
#seqtab.roots, seqtab.soil
seqtab.roots <- makeSequenceTable(mergers)
seqtab.soil <- makeSequenceTable(mergers)

#Track reads through pipeline
getN <- function(x) sum(getUniques(x))
trackroots2875<-cbind(outroots2875, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
trackroots5675<-cbind(outroots5675, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
tracksoil5180<-cbind(outsoil5180, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) 
head(trackroots2875)
head(trackroots5675)
head(tracksoil5180)
hist(trackroots5675[,5]-trackroots2875[,5],breaks=65)
sort(trackroots4447[,5]-trackroots7570[,5])


###### MERGE 2 sets ######
seqtab.all<-mergeSequenceTables(seqtab.roots,seqtab.soil)

#remove chimeras, start 11:30, end 11:31
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim5675 <- removeBimeraDenovo(seqtab.roots, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim2875
seqtab.nochim5180s <- removeBimeraDenovo(seqtab.soil, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim5675[1:5,1:5]
seqtab.nochim2875[1:5,1:5]


write.csv(seqtab.nochim, "seqtab.nochim.csv")
rownames(seqtab.nochim)

#inspect distribution of sequence length
plot(table(nchar(getSequences(seqtab.nochim5180s))))
table(nchar(getSequences(seqtab.nochim5180s)))
#Roots
#using the first 2875 trimming it was 20960 are 254, 544 are 255, max 460 (but only 1 sequence)
#using 5675 trimming it was 19799 are 253, 2157 are 254, max 440 (but only 1 sequence). When I looked at the sequences the 2875 ones were 254 and the last base pair is an A, the 5675 ones were 253 and were identical to the bp 1-253 from the first try. Thus it is true that the first run is just adding an A to the end of the sequences, this A is from the reverse primer. OOH It it because the reads are actually 301 bp, not 300!!!

#Soils:
#using 5180 trimming it was 36121 are 253, 3832 are 254, max 457 (only 2 sequences)

#normally do track here but not if you merged multiple lanes
#trackroots <- cbind(outroots2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))




###### Assign taxonomy ######

# to get the fasta file gohere:https://unite.ut.ee/repository.php and be sure to do "general release fasta" not qiime. This is 4/4/24 releae, the same one I used for the mangrove project. Note that the zip files is named slightly differently than the actual fasta when you unzip it (_s_ vs _dyanmic_s_)
#unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Survey/Stats/Gradient/QIIME2/sh_general_release_dynamic_s_04.02.2020.fasta" #
#unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/Unite/sh_general_release_dynamic_25.07.2023.fasta" #
unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Mangroves/Stats/sh_general_release_dynamic_s_04.04.2024.fasta"

#start 2:23, done  (70 min)
taxa<-assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, minBoot=70, tryRC = TRUE,outputBootstraps=T)
taxaonly<-taxa$tax

write.csv(taxa,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/taxa.csv")


###### Create tax table for Phyloseq ######

sequences<-as.data.frame(rownames(taxaonly))
rownames(sequences)<-paste0("OTU",1:nrow(sequences))
sequences$OTU<-rownames(sequences)
colnames(sequences)<-c("sequence","OTUID")
#OTUID<-as.data.frame(sequences$OTUID)

taxa1<-cbind(as.data.frame(sequences$OTUID), taxaonly)
write.csv(taxa1,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/taxa_OTUID.csv")

#OTUID<-sequences$OTUID
rownames(taxa1)<-sequences$OTUID

#use this for phyloseq tax table
taxa1<-taxa1[,2:8]

phy.tax<-tax_table(as.matrix(taxa1))
phy.tax



###### Create sample data table ######
sam<-as.data.frame(rownames(seqtab.nochim))
names(sam)<-"SampleNumber"

sam2<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Data/NWTrootssoil2023.csv") #roots only is NWTroots2023.csv

#join
sam3<-left_join(sam,sam2, by="SampleNumber")

#add column for decontam and clean
sam3$Sample_or_Control<-ifelse(sam3$PlotID=="control","Control","Sample")
sam3$is.neg<-ifelse(sam3$Sample_or_Control=="Control",TRUE,FALSE)

#start to make it a phyloseq object
phy.sam<-sample_data(sam3)
rownames(phy.sam)<-sam3$SampleNumber


###### Create official ASV Table ######
OTU.Table<-as.data.frame(seqtab.nochim)
colnames(OTU.Table)<-sequences$OTUID
#temp<-1:dim(OTU.Table)[2] #16535 OTUs
#colnames(OTU.Table)<-paste("OTU",temp,sep="")
OTU.Table[1:5,1:5]


#check that otu table is in the same order as sample file
cbind(rownames(OTU.Table),as.character(sam3$SampleNumber))
#rownames(OTU.Table)<-SampleNumber
#OTU.Table$SampleID_Original<-rownames(OTU.Table)
#ASV.Table<-left_join(OTU.Table,sam,by="SampleID_Original")
#rownames(ASV.Table)<-ASV.Table$SampleID
#ASV.Table<-ASV.Table[,1:16456]

phy.ASV<-otu_table(OTU.Table,taxa_are_rows = F)
phy.ASV[,1:5]

#create phyloseq object
phy<-merge_phyloseq(phy.ASV,phy.sam,phy.tax)
sample_data(phy)
otu_table(phy)[1:5,1:6]


#### Decontam ####

datITS<-phy

print(as.data.frame(unique(tax_table(datITS)[,"Kingdom"])),row.names=F)

####Filter and remove contaminant sequences with decontam()####
#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement

#Filter samples with really low numbers of reads (<1000), not sure if this is necessary for the prevalence method but it is kind of weird that some of my actual samples had lower number of reads than the negative controls. This could be because the pcr failed and we put them in anyway "just in case" but if the pcr failed and all we got was contamination, then that skews the ability to see that they are contaminated. here it was only 1 sample
sort(sample_sums(datITS))

datITSS<-datITS%>%
  subset_samples(sample_sums(datITS)>900) 
datITSS

sample_data(datITSS)

dfITS <- as.data.frame(sample_data(datITSS)) # Put sample_data into a ggplot-friendly data.frame
dfITS$LibrarySize <- sample_sums(datITSS)
dfITS <- dfITS[order(dfITS$LibrarySize),]
dfITS$Index <- seq(nrow(dfITS))
ggplot(data=dfITS, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

contamdf.prevITS <- isContaminant(datITSS, method="prevalence", neg="is.neg")
table(contamdf.prevITS$contaminant)#default threshold is .1
#contamdf.prevITS5 <- isContaminant(datITSS, method="prevalence", neg="is.neg",threshold = .2)
#table(contamdf.prevITS5$contaminant)

datITSS.pa <- transform_sample_counts(datITSS, function(abund) 1*(abund>0))
datITSS.pa.neg <- prune_samples(sample_data(datITSS.pa)$is.neg == TRUE, datITSS.pa)
datITSS.pa.pos <- prune_samples(sample_data(datITSS.pa)$is.neg == FALSE, datITSS.pa)

# Make data.frame of prevalence in positive and negative samples
#threshold .1
datITSSdf.pa <- data.frame(pa.pos=taxa_sums(datITSS.pa.pos), pa.neg=taxa_sums(datITSS.pa.neg),contaminant=contamdf.prevITS$contaminant)
ggplot(data=datITSSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#threshold .2
#datITSSdf.pa <- data.frame(pa.pos=taxa_sums(datITSS.pa.pos), pa.neg=taxa_sums(datITSS.pa.neg),contaminant=contamdf.prevITS5$contaminant)
#ggplot(data=datITSSdf.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#I will use a threshold = 0.1 (in phrag suff I think I've used .2). (when I did roots alone I used 0.1) 

#take out contaminants and filter negative controls, there were 98 contaminants
datITSS2 <- prune_taxa(!contamdf.prevITS$contaminant, datITSS)
datITSS3 <-datITSS2 %>%
  subset_samples(is.neg==FALSE)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T) #odd some taxa additionally are removed after removing the contaminants. these must be taxa only found in the negative controls (which for some reason are not called contaminants) (??)

min(sample_sums(datITSS4))
sort(sample_sums(datITSS3))

#Filter out root samples and samples with low sampling depth. THis is tough, I could rarefy to 4476 and only delete 4 low samples or I could rarefy to 6131 and delete 6 low samples. I will choose 4476 so that I can have more root and soil samples from the same plots
datITSS4 <-datITSS3 %>%
  #subset_samples(SampleType=="soil")%>%
  subset_samples(sample_sums(datITSS3)>4000) %>%
  filter_taxa(function(x) sum(x) > (0), prune=T)


####### Rarefy to 4476 ######
#This takes out samples: r138, r137, r139 (this was taken out before decontam), s23,  s48
datITSS4
datITSS5<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
datITSS5c<-datITSS4%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS4)),rngseed=10,replace=F)
#2320 OTUs were removed because they are no longer present in any sample after random subsampling




#Adding to the sample data with dat from dat17
#need to open and run LoadingPlantSoilData.R script

# tempf<-sample_data(datITSS5c)
# tempf$Plot<-as.numeric(as.character(tempf$Plot))
# tempd<-dat17[,c(1:4,52:77)]
# 
# tempb<-left_join(tempf,tempd,by="Plot")#warning is ok
# tempb2<-cbind(Plottest=tempb$Plot,tempb[,6:34])
# 
# sample_data(datITSS5c)[,6:35]<-tempb2
# sample_data(datITSS5c)[1:5,1:35]
# 
# #checks
# cbind(as.character(sample_data(datITSS5c)$SampleID),as.numeric(as.character(sample_data(datITSS5c)$Plot)),sample_data(datITSS5c)$Plottest)
# sample_data(datITSS5c)$Plottest<-NULL
# 
# sample_data(datITSS5)
# sample_data(datITSS5)[,6:35]<-tempb2
# cbind(as.numeric(as.character(sample_data(datITSS5)$Plot)),sample_data(datITSS5)$Plottest)
# sample_data(datITSS5)$Plottest<-NULL


#Make OTU tables, this takes a while, start 10:47-10:52
datITSS5otu<-cbind(sample_data(datITSS5),otu_table(datITSS5))
datITSS5cotu<-cbind(sample_data(datITSS5c),otu_table(datITSS5c))




#Richness

richITS<-estimate_richness(datITSS5c, split = TRUE, measures = c("Observed", "Shannon","Chao1","Simpson","InvSimpson"))
colnames(richITS)<-c("RichnessITS","Chao1ITS","se.chao1ITS","ShannonITS","SimpsonITS","InvSimpsonITS")
richITS$SampleID<-rownames(richITS)
#richITS2<-separate(richITS,SampleID,c(NA,"Plot"),"s")
#richITS2$Plot<-as.numeric(richITS2$Plot)
#richITS2<-richITS%>%
#  left_join(data.frame(sample_data(datITS2rcsoil)))
#sample_data(datITSS5c)[,]

dat17b<-dat17%>%
  full_join(richITS2)
dat17<-dat17b




###### FUNguild ######
#https://github.com/UMNFuN/FUNGuild
#https://github.com/brendanf/FUNGuildR
devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)

#The weird thing about FUNGuildR is that for some taxa, like when I have Fusarium equiseti it matches to Nectriaceae rather than Fusarium in the database. When I used the FUNGUILD in the terminal, it worked better and matched to the lowest taxonomic level. 

#Trying with FUNguildR
# funguild<-data.frame(tax_table(datITSS5c))
# funguild2<-funguild%>%
#   unite("Kingdom_Species",sep=";",remove=T)
# colnames(funguild2)<-"Taxonomy"
# funguild2$OTUID<-rownames(funguild2)
# funguild3<-data.frame(OTUID=funguild2$OTUID,Taxonomy=funguild2$Taxonomy)
# head(funguild3)
# funguild4<-funguild_assign(funguild3)
# View(funguild4)

#I don't need this but I'm curious if it differs from the online one here: http://www.stbates.org/funguild_db.php   it doesn't seem to differ...
# fung <- get_funguild_db()
# fung[fung$taxon=="Meliniomyces",]

#Trying in terminal
temp<-data.frame(tax_table(datITSS5c))
temp2<-temp%>%
  unite("Kingdom_Species",sep=";",remove=T)
head(temp2)
colnames(temp2)<-"taxonomy"

datITSS5cfunguild<-cbind(t(otu_table(datITSS5c)),temp2)
datITSS5cfunguild[1:5,155:157]
datITSS5cfunguild[1:5,1:5]
write.csv(datITSS5cfunguild,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/datITSS5cfunguild.csv")
#then open the file and write "OTU ID" in the first cell

#then open terminal and navigate to folder and run
#note when you download the py file from github make sure you don't just download (save as) the html file, make sure it's a python file, open the code and click the download raw file button
python Guilds_v1.1.py -otu datITSS5cfunguild.csv -m -u

#then open datITSS5cfunguild.guilds.txt and delete OTU ID in the first cell

funguild<-read.delim2("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot/NiwotIndirectEffects/Stats/QIIME2/datITSS5cfunguild.guilds.txt",header=T,row.names = 1)
head(funguild)
funguild[20:25,300:323]
funguild$OTU<-rownames(funguild)
funguild
funguild$OTUnum <- as.numeric(gsub("^.{0,3}", "", funguild$OTU))
funguild2<-funguild%>%
  arrange(OTUnum)
funguild2[1:5,320:325]
funguild2[1:5,1:5]

#comparing the two, often the results for guild or trophic mode aren't different, but sometimes the FUNguildR result is actually more specific/narrowed down which is weird, but I guess the database does not automatically have the family have all the charactertics of all genera. so if not a lot is know just about a family you might get only one trophic mode, whereas for a genus you get two b/c there are more studies at the genus level. this is odd behavior but I think doing it to genus is better
# cbind(funguild2$Taxon,funguild4$taxon)
# i<-258
# funguild2$taxonomy[i]
# funguild4$Taxonomy[i]
# funguild2$Taxon[i]
# funguild4$taxon[i]
# funguild2$Guild[i] #terminal
# funguild4$guild[i] #FUNguildR
# funguild2$Trophic.Mode[i] #terminal
# funguild4$trophicMode[i] #FUNguildR
# funguild2[i,]
# funguild4[i,]




#rarefied to 4476
colSums(funguild2[,2:5])
funguild2[c("Trophic.Mode","Confidence.Ranking")]
funguild2[1:5,315:325]

funguild2%>%filter(Trophic.Mode==("Pathotroph"))
funguild2%>%filter(Trophic.Mode==("Symbiotroph"))

#funguild_query doesn't work on the terminal output b/c the column names are slightly different from the FUNguildR output
#funguild_query("*Saprotroph*", "Trophic.Mode", db = funguild2)

funguild3<-funguild2%>%
  filter(Trophic.Mode%in%c("Symbiotroph","Pathotroph","Saprotroph"))%>%
  #filter(Confidence.Ranking%in%c("Probable","Highly Probable"))%>%
  arrange(Trophic.Mode)%>%
  dplyr::select(r1:s99,Trophic.Mode)%>%
  group_by(Trophic.Mode)%>%
  summarise_all(list(sum=sum))
funguild4<-data.frame(funguild3)
row.names(funguild4)<-funguild3$Trophic.Mode;funguild4$Trophic.Mode<-NULL
funguild5<-data.frame(t(funguild4))
head(funguild5)
funguild5$SampleNumbertemp<-row.names(funguild5)
funguild6<-funguild5%>%
  separate(SampleNumbertemp,into=c("SampleNumber",NA))
#funguild6$SampleNumber<-as.numeric(sub("X","",funguild6$SampleNumber))
funguild6[,1:3]<-funguild6[,1:3]/4476*100

temp<-sample_data(datITSS5c)
cbind(temp$SampleNumber,funguild5$SampleNumber)
temp[,18:20]<-funguild6[,1:3]
temp<-data.frame(temp)

#Make a graph for annual report
#pathogens
m1<-temp%>%
  filter(PlotType!="Survey")%>%
  group_by(SampleType,CommunityType,Treatment)%>%
  summarise(mean=mean(Pathotroph), se=std.error(Pathotroph))
m1
m1$CommunityType<-factor(m1$CommunityType,levels=c("WM","MM","DM"))
#m1$CommunityType<-recode_factor(m1$CommunityType,"C"="Control","E"="Treatment")

pdf("Figs/PathogensExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Treatment, y=mean,color=SampleType))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Pathogens %")+
  #  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(CommunityType),strip.position = "bottom")
dev.off()

#symbionts
m1<-temp%>%
  filter(PlotType!="Survey")%>%
  group_by(SampleType,CommunityType,Treatment)%>%
  summarise(mean=mean(Symbiotroph), se=std.error(Symbiotroph))
m1
m1$CommunityType<-factor(m1$CommunityType,levels=c("WM","MM","DM"))
#m1$CommunityType<-recode_factor(m1$CommunityType,"C"="Control","E"="Treatment")

pdf("Figs/SymbiontsExperiment.pdf",width=3.2,height=2.2)
ggplot(data=m1, aes(x=Treatment, y=mean,color=SampleType))+   
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.2,size=.5)+
  geom_point(size=1.8,show.legend = FALSE)+#, aes(group=Seed.Origin, fill=Seed.Origin, shape=Seed.Origin, color = Seed.Origin)
  ylab("Symbionts %")+
  #  ylim(5,60)+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=9),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside",axis.title.x = element_blank())+
  facet_wrap(vars(CommunityType),strip.position = "bottom")
dev.off()



dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#To get "Plant Pathogen" or "Endophyte" or "Arbuscular Mycorrhizal"
funguild2<-funguild%>%
  filter(Guild%in%c("Plant Pathogen","Endophyte","Arbuscular Mycorrhizal"))%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1:3]<-funguild5[,1:3]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#Note: any taxon classified as one thing (i.e. plant pathogen with nothing else) always gets at least a probable, no possibles.
#Trying taking out the "probables" for plant pathogens so it only includes highly probable. if i do this, there are only 7 non zero plots for plant pathogens, so that's not reasonable 
funguild2<-funguild%>% 
  filter(Guild%in%c("Plant Pathogen","Endophyte","Arbuscular Mycorrhizal"))%>%
  filter(Confidence.Ranking=="Highly Probable")%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1:3]<-funguild5[,1:3]/5984*100
colnames(funguild5)[1:3]<-c("Arbuscular.Mycorrhizalhighlyprobable","Endophytehighlyprobable","Plant.Pathogenhighlyprobable")

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

#Getting data for ectomycorrhizal
funguild2<-funguild%>% 
  filter(Guild%in%c("Ectomycorrhizal"))%>%
  arrange(Guild)%>%
  dplyr::select(s1:s98,Guild)%>%
  group_by(Guild)%>%
  summarise_all(list(sum=sum))
funguild3<-data.frame(funguild2)
row.names(funguild3)<-funguild3$Guild;funguild3$Guild<-NULL
funguild4<-data.frame(t(funguild3))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1]<-funguild5[,1]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b

dat17$PlantMutualist<-dat17$Ectomycorrhizal+dat17$Arbuscular.Mycorrhizal+dat17$Endophyte


#Try any description that includes "plant pathogen"
ind<-grep("Plant Pathogen",funguild$Guild)
funguild2<-funguild[ind,]
funguild3<-colSums(funguild2[,1:162])
funguild4<-data.frame(t(data.frame(t(funguild3))))
head(funguild4)
funguild4$Plottemp<-row.names(funguild4)
colnames(funguild4)[1]<-"Plant.Pathogenanywhere"
funguild5<-funguild4%>%
  separate(Plottemp,into=c("Plot",NA))
funguild5$Plot<-as.numeric(sub("s","",funguild5$Plot))
funguild5[,1]<-funguild5[,1]/5984*100

dat17b<-dat17%>%
  full_join(funguild5)
dat17<-dat17b
