#Network analysis

##### Fungi #####
tempnet<-datITSS5cE%>%
  subset_samples(SampleType=="soil")%>%
  #subset_samples(Site=="EastKnoll")%>%
  #subset_samples(Treatment=="Experimental")%>%
  subset_samples(Treatment=="Control")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)%>%
  filter_taxa(function(x) sum(x) > 200, prune=T) #40
tempnet

#When doing a single site (low sample size), I think it is best to use spearman and a t-test, otherwise you get a million edges
net_pears <- netConstruct(tempnet,  
                          measure = "spearman",#pearson
                          normMethod = "clr",
                          zeroMethod = "pseudo",#multRepl #pseudo just adds 1 to all counts
                          sparsMethod = "t-test",#threshold
                          alpha=0.10,
                          thresh = 0.3,
                          verbose = 3)

#When doing it on all sites, pearson seems better, with alpha=0.05
net_pears <- netConstruct(tempnet,  
                          measure = "pearson",#
                          normMethod = "clr",
                          zeroMethod = "pseudo",#multRepl #pseudo just adds 1 to all counts
                          sparsMethod = "t-test",#threshold
                          alpha=0.05,
                          thresh = 0.3,
                          verbose = 3)

#Takes 10 minutes with 500 taxa
props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")

#With stricter edge filter
plot(props_pears,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.4,
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     #labelScale = FALSE,
     #cexLabels = 1.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = paste0("Network on OTU level with Pearson correlations",
                     "\n(edge filter: threshold = 0.4)"),
     showTitle = F,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)

plot(props_pears, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

#Prettier
plot(props_pears, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     #labelScale = FALSE,
     #cexLabels = 1.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)




###### Network comparison ######

tempnet<-datITSS5cE%>%
  subset_samples(SampleType=="soil")


# Split the phyloseq object into two groups
tempnet_co <- phyloseq::subset_samples(tempnet, Treatment == "Control")
tempnet_ex <- phyloseq::subset_samples(tempnet, Treatment == "Experimental")

tempnet_co
tempnet_ex

#set sample size of smaller group
n_yes <- phyloseq::nsamples(tempnet_ex)

# Network construction


#from createassocperm
net_amgut <- netConstruct(counts_matched, 
                          group = group_vec,
                          matchDesign = c(1,2),
                          filtTax = "highestFreq",
                          filtTaxPar = list(highestFreq = 50),
                          measure = "pearson",
                          zeroMethod = "pseudo", 
                          normMethod = "clr",
                          sparsMethod = "threshold", 
                          thresh = 0.4,
                          seed = 123456)

# Network analysis with default values
props_amgut <- netAnalyze(net_amgut)






###### From Susannah's code ######
net_mh <- netConstruct(data = tempnet_co, 
                       data2 = tempnet_ex,
                       filtTax = "highestFreq",
                       filtTaxPar = list(highestFreq = 500), # cite dini-andreote
                       zeroMethod = "pseudo", 
                       normMethod = "clr", #center log-ratio transformation
                       measure = "sparcc", #pearson
                       sparsMethod = "threshold", #t-test
                       thresh = 0.4, #threshold for sparsification
                       cores = 5,
                       seed = 2024)

#500 makes things really slow
net_mhj <- netConstruct(data = tempnet_co, 
                       data2 = tempnet_ex,
                       jointPrepro = T, #is the filtering done on the whole dataset or each set individually
                       filtTax = "highestFreq",
                       filtTaxPar = list(highestFreq = 500), # cite dini-andreote
                       zeroMethod = "pseudo", 
                       normMethod = "clr", #center log-ratio transformation
                       measure = "sparcc", #pearson
                       sparsMethod = "threshold", #t-test
                       thresh = 0.4, #threshold for sparsification
                       cores = 5,
                       seed = 2024)

props_mh <- netAnalyze(net_mh,clustMethod = "cluster_fast_greedy")
props_mhj <- netAnalyze(net_mhj,clustMethod = "cluster_fast_greedy")

#?plot.microNetProps
plot(props_mhj, 
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.3,
     sameLayout = T, #can change to F
     repulsion = 0.9,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "degree", #mclr
     nodeFilter = "clustMin",
     nodeFilterPar = 10,
     nodeSizeSpread = 3, 
     edgeTranspLow = 50,
     edgeTranspHigh = 30,
     labels=F,
     labelScale = FALSE,
     nodeTransp = 30,
     cexNodes = 1.5, 
     cexLabels = 0.7,
     cexHubLabels = 1,
     cexTitle = 1,
     groupNames = c("Control", "Experimental"),
     hubBorderCol  = "gray40")

#start 2:00 - 2:09, I have 8 cores (6 performance, 2 efficiency)
comp_mh <- netCompare(props_mh, permTest = TRUE, nPerm = 20,#nPerm=1000 
                      storeAssoPerm = TRUE,
                      fileStoreAssoPerm = "assoPerm_comp_mh",
                      storeCountsPerm = FALSE, 
                      testRand = F,#was true
                      cores = 5,
                      seed = 123456)
#start 8:20pm, 80% at 11:20. it broke connection and stopped at 80% maybe when I saved the script file (?)
comp_mhj <- netCompare(props_mhj, permTest = TRUE, nPerm = 100,#nPerm=1000 
                      storeAssoPerm = TRUE,
                      fileStoreAssoPerm = "assoPerm_comp_mhj",
                      storeCountsPerm = FALSE, 
                      testRand = F,#was true
                      cores = 5,
                      seed = 123456)
summary(comp_mh, groupNames=c("Control", "Experimental"))
summary(comp_mhj, groupNames=c("Control", "Experimental"))






#From net_comparison

net_season <- netConstruct(data = tempnet_co, 
                           data2 = tempnet_ex,  
                           filtTax = "highestVar",
                           filtTaxPar = list(highestFreq = 500),#highestVar = 50
                           filtSamp = "highestFreq",
                           filtSampPar = list(highestFreq = n_yes),
                           measure = "spring",
                           measurePar = list(nlambda = 10, 
                                             rep.num = 10,
                                             Rmethod = "approx"),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 123456)

props_season <- netAnalyze(net_season, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)
summary(props_season)

plot(props_season, 
     sameLayout = TRUE, 
     repulsion = 0.95,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     labelScale = FALSE,
     cexNodes = 1.5, 
     cexLabels = 1,
     cexHubLabels = 1,
     cexTitle = 1,
     groupNames = c("Control", "Experimental"),
     hubBorderCol  = "gray40")







##### Bacteria #####
tempnetb<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")%>%
  #subset_samples(Site=="EastKnoll")%>%
  subset_samples(Treatment=="Experimental")%>%
  #subset_samples(Treatment=="Control")%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)%>%
  filter_taxa(function(x) sum(x) > 100, prune=T) #40
tempnetb #shoot for like 100-150 taxa


#When doing it on all sites, pearson seems better, with alpha=0.05
net_pearsb <- netConstruct(tempnetb,  
                          measure = "pearson",#
                          normMethod = "clr",
                          zeroMethod = "pseudo",#multRepl #pseudo just adds 1 to all counts
                          sparsMethod = "t-test",#threshold
                          alpha=0.05,
                          thresh = 0.3,
                          verbose = 3)

#Takes 10 minutes with 500 taxa
props_pearsb <- netAnalyze(net_pearsb, clustMethod = "cluster_fast_greedy")

#With stricter edge filter
plot(props_pearsb,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.6,
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.7,
     rmSingles = TRUE,
     #labelScale = FALSE,
     #cexLabels = 1.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = paste0(""),
     showTitle = F,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)





###### Network comparison ######

tempnetb<-dat16SS5cE%>%
  subset_samples(SampleType=="soil")
tempnetb

# Split the phyloseq object into two groups
tempnetb_co <- phyloseq::subset_samples(tempnetb, Treatment == "Control")
tempnetb_ex <- phyloseq::subset_samples(tempnetb, Treatment == "Experimental")

tempnetb_co
tempnetb_ex

#set sample size of smaller group
n_yes <- phyloseq::nsamples(tempnetb_ex)



# Network construction

###### From Susannah's code ######
net_mhb <- netConstruct(data = tempnetb_co, 
                       data2 = tempnetb_ex,
                       filtTax = "highestFreq",
                       filtTaxPar = list(highestFreq = 500), # cite dini-andreote
                       zeroMethod = "pseudo", 
                       normMethod = "clr", #center log-ratio transformation
                       measure = "sparcc", #pearson
                       sparsMethod = "threshold", #t-test
                       thresh = 0.4, #threshold for sparsification
                       cores = 5,
                       seed = 2024)

#500 makes things really slow
net_mhjb <- netConstruct(data = tempnetb_co, 
                        data2 = tempnetb_ex,
                        jointPrepro = T, #is the filtering done on the whole dataset or each set individually
                        filtTax = "highestFreq",
                        filtTaxPar = list(highestFreq = 200), # cite dini-andreote
                        zeroMethod = "pseudo", 
                        normMethod = "clr", #center log-ratio transformation
                        measure = "sparcc", #pearson
                        sparsMethod = "threshold", #t-test
                        thresh = 0.4, #threshold for sparsification
                        cores = 5,
                        seed = 2024)

props_mhb <- netAnalyze(net_mhb,clustMethod = "cluster_fast_greedy")
props_mhjb <- netAnalyze(net_mhjb,clustMethod = "cluster_fast_greedy")

#?plot.microNetProps
plot(props_mhjb, 
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.3,
     sameLayout = T, #can change to F
     repulsion = 0.9,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "degree", #mclr
     nodeFilter = "clustMin",
     nodeFilterPar = 10,
     nodeSizeSpread = 3, 
     edgeTranspLow = 50,
     edgeTranspHigh = 30,
     labels=F,
     labelScale = FALSE,
     nodeTransp = 30,
     cexNodes = 1.5, 
     cexLabels = 0.7,
     cexHubLabels = 1,
     cexTitle = 1,
     groupNames = c("Control", "Experimental"),
     hubBorderCol  = "gray40")

#start 2:00 - 2:09, I have 8 cores (6 performance, 2 efficiency)
comp_mhb <- netCompare(props_mhb, permTest = TRUE, nPerm = 20,#nPerm=1000 
                      storeAssoPerm = TRUE,
                      fileStoreAssoPerm = "assoPerm_comp_mhb",
                      storeCountsPerm = FALSE, 
                      testRand = F,#was true
                      cores = 5,
                      seed = 123456)
#start 11:27, end 12:20 (200 taxa)
comp_mhjb <- netCompare(props_mhjb, permTest = TRUE, nPerm = 100,#nPerm=1000 
                       storeAssoPerm = TRUE,
                       fileStoreAssoPerm = "assoPerm_comp_mhjb",
                       storeCountsPerm = FALSE, 
                       testRand = T,#was true
                       cores = 5,
                       seed = 123456)
summary(comp_mhb, groupNames=c("Control", "Experimental"))
summary(comp_mhjb, groupNames=c("Control", "Experimental"))



