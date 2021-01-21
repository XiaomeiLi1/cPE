## Data dived PPI but filted by  the directed PPI

# This is the script of the proposed method and its application in detecting Eo and Lo preeclampsia biomarkers. 
# The method is applied to the preeclampsia dataset to discover Eo and Lo preeclampsia biomarkers.
# 
# To run the script, please prepare below input files and reset environment variables in the script.
# 
# The input files include the followings and are put in the folder "rootDir/Data" (rootDir is an environment variable):
#   
# * PPI.xls - Protein protein interaction network
# 
# * Browse Transcription Factors hg19 - resource_browser.csv - Transcription factors (TFs)
# 
# * GSE75010.rda - PE expression data
# 
# 
# This script uses a library from the paper Liu, Y.-Y., et al. (2011). "Controllability of complex networks." Nature 473: 167. The code of the library can be downloaded from https://scholar.harvard.edu/yyl/code. You need to build the code before running this script.

#=========================================================================
#=========================================================================
# Eo and Lo preeclampsia biomarkers
#=========================================================================
#=========================================================================
# Clear the environment if needed
#rm(list = ls())

# Load necessary libraries if any
library(readxl)
library(ggplot2)
library(varhandle)
library(scales)
library(reshape)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(xtable)

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
rootDir <- "C:/unisa/RstudioProjects/PEDriver" # And put the input files in "rootDir/Data"
controlDir <- "C:/unisa/RstudioProjects/PEDriver/External" # Put here the library from 
# the paper "Controllability of complex networks."

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
# # For EoPE biomarkers
outDir <- "C:/unisa/RstudioProjects/PEDriver/Data/Output/EoPE" # Output folder
type <- "EoPE"
#---------------------------------------

# Include the script of functions
source(paste(rootDir, "/R/ProposedMethod_Functions.R", sep=""))

# Main script 
#================================================================
# (1) Building the network for a specific condition
#================================================================
# Load the tumor expression data
#load(paste(rootDir, "/Data/GSE75010.rda", sep = ""))

# the data is a list(length = 2), contain PEType, mRNAs
# $mRNAs is sample*mRNA with rownames and colnames
PE_matchedData = list(length = 2)
PE_matchedData$PEType = type
index = which(pd$subtype == type)
ddata = exprs(GSE75010)
PE_matchedData$mRNAs = t(ddata[,index])

# Get PPI network to a dataframe (34814*5)
edges <- read_excel(paste(rootDir, "/Data/PPI.xls",
                          sep = ""), sheet = 1)
interactions <- edges[, c(1, 3)] # gene symbol
colnames(interactions) <- c("cause", "effect")
# which interactions exist in the current dataset
interactions <- interactions[which(interactions$cause %in% colnames(PE_matchedData$mRNAs)),]
interactions <- interactions[which(interactions$effect %in% colnames(PE_matchedData$mRNAs)),]
nodes <- unique(union(interactions$cause, interactions$effect))

# TFs: Download the TF list from http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19
tfs <- read.csv(paste(rootDir, "/Data/Browse Transcription Factors hg19 - resource_browser.csv",
                      sep = ""))
i <- which(levels(tfs$Symbol) %in% nodes)
# extract TF exprssion data
tfData <- PE_matchedData$mRNAs[, levels(tfs$Symbol)[i]]

# Update PE data of mRNAs, only contain mRNAs in the network and not in the TF list
PE_matchedData$mRNAs <- PE_matchedData$mRNAs[,
                                             nodes[which(!(nodes %in% levels(tfs$Symbol)[i]))]]
mRNAsData_PE <-  PE_matchedData$mRNAs

# Combine data to a single matrix with columns are miRNA, mRNA, TFs
nomR <- ncol(PE_matchedData$mRNAs)
noTF <- ncol(tfData)
PE_data <- cbind(mRNAsData_PE, tfData)

# Free the memory
gc()

# Build the network
PE_network <-  buildNetworkForPersonalised(interactions, nomR, noTF, PE_data, 
                                           rootDir,usingpVal = TRUE, cutoff = 0.05)

# Save the network
write.csv(PE_network, paste(outDir, "/pVal_PE_network.csv", sep = ""), row.names = FALSE)

# Analyse network summary of the network
# PE_network <- read.csv(paste(outDir, "/PE_network.csv", sep = ""))

analyseNetworkForPersonalised(nomR, noTF, PE_network, PE_data, 
                              paste(outDir, "/pVal_PE_network_analysis.txt", sep = ""))

#================================================================
# (2) Identifying PE drivers
#================================================================

# Analyse controllability of the network
interactions <- read.csv(paste(outDir, "/pVal_PE_network.csv",
                               sep = ""))
# Write the edges of the network for analysing controllability
write.table(interactions, paste(outDir, "/Controllability/pVal_edges.dat", sep = ""), row.names = FALSE, col.names=FALSE, quote=FALSE)
# Run the controllability analysis
cmd <- paste(controlDir, "/parse.exe ", outDir, 
             "/Controllability/pVal_edges.dat", sep = "")
system(cmd)
cmd <- paste(controlDir, "/ControllabilityAnalysis.exe ", outDir, 
             "/Controllability/pVal_edges.dat", sep = "")
system(cmd)
# Analyse controllability of the network and output in a file
analyseControllability(paste(outDir, "/Controllability/pVal_edges.dat.output", sep = ""),
                       paste(outDir, "/pVal_analyseControllability.txt", sep = ""))

# Identify critical nodes in the network
# Read the result
nodetype <- read.table(paste(outDir, "/Controllability/pVal_edges.dat.nodetype", sep = ""))
colnames(nodetype) <- c("Name", "K", "Kin", "Kout", "TypeI", "TypeII")
nodetype$Name = as.character(nodetype$Name)
# Type 1 critical nodes of the network
Type1_critical_nodes <- nodetype[which(nodetype$TypeII == 0),]
Type1_critical_nodes$Name =  as.character(Type1_critical_nodes$Name)
Type2_critical_nodes <- nodetype[which(nodetype$TypeI == 0),]
Type2_critical_nodes$Name =  as.character(Type2_critical_nodes$Name)
# Save file
write.csv(Type1_critical_nodes, paste(outDir, "/pVal_Type1_critical_nodes.csv", sep = ""),
          row.names = FALSE)
write.csv(Type2_critical_nodes, paste(outDir, "/pVal_Type2_critical_nodes.csv", sep = ""),
          row.names = FALSE)

# Combine Control data, rank candidate PE drivers
index = c(which(pd$subtype == type),which(pd$subtype == "Control" & pd$GA < 34))

index = c(which(pd$subtype == type),which(pd$phenotype == "normal" & pd$GA >= 34))

temp = pd$subtype[index]
temp[temp=="NA"] <- "Control"
design <- model.matrix(~0+factor(temp))
colnames(design)=levels(factor(temp))
rownames(design)=colnames(ddata[,index])
design

contrast.matrix<-makeContrasts(paste0(unique(temp),collapse = "-"),levels = design)
contrast.matrix 

##step1
fit <- lmFit(ddata[,index],design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
limma_notrend_results_EoPE = na.omit(tempOutput) 
head(limma_notrend_results_EoPE)
write.csv(limma_notrend_results_EoPE, paste0(outDir, "/limma_notrend_results_EoPE.csv"),quote = F)

limma_notrend_results_EoPE = read.csv(paste0(outDir, "/limma_notrend_results_EoPE.csv"),header = T, stringsAsFactors = F)
colnames(limma_notrend_results_EoPE)[1] = "gene"
rownames(limma_notrend_results_EoPE) = limma_notrend_results_EoPE$gene

limma_notrend_results_EoPE$gene = rownames(limma_notrend_results_EoPE)
limma_notrend_results_EoPE = limma_notrend_results[which(limma_notrend_results_EoPE$adj.P.Val < 0.05),c(1,2)]
limma_notrend_results_EoPE$logFC = abs(limma_notrend_results_EoPE$logFC)

Type1_candidate_EoPE_drivers <- merge(Type1_critical_nodes, limma_notrend_results_EoPE,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")

Type2_candidate_EoPE_drivers <- merge(Type2_critical_nodes, limma_notrend_results,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")


# By the absolute log-fold-change
Type1_candidate_EoPE_drivers <- 
  Type1_candidate_EoPE_drivers[order(Type1_candidate_EoPE_drivers$logFC, decreasing = TRUE),]

Type2_candidate_EoPE_drivers <- 
  Type2_candidate_EoPE_drivers[order(Type2_candidate_EoPE_drivers$logFC, decreasing = TRUE),]

# # For LoPE biomarkers
outDir <- "C:/unisa/RstudioProjects/PEDriver/Data/Output/LoPE" # Output folder
type <- "LoPE"
#---------------------------------------

# Include the script of functions
# source(paste(rootDir, "/R/ProposedMethod_Functions.R", sep=""))

# Main script 
#================================================================
# (1) Building the network for a specific condition
#================================================================
# Load the tumor expression data
#load(paste(rootDir, "/Data/GSE75010.rda", sep = ""))

# the data is a list(length = 2), contain PEType, mRNAs
# $mRNAs is sample*mRNA with rownames and colnames
PE_matchedData = list(length = 2)
PE_matchedData$PEType = type
index = which(pd$subtype == type)
ddata = exprs(GSE75010)
PE_matchedData$mRNAs = t(ddata[,index])

# Get PPI network to a dataframe (34814*5)
edges <- read_excel(paste(rootDir, "/Data/PPI.xls",
                          sep = ""), sheet = 1)
interactions <- edges[, c(1, 3)] # gene symbol
colnames(interactions) <- c("cause", "effect")
# which interactions exist in the current dataset
interactions <- interactions[which(interactions$cause %in% colnames(PE_matchedData$mRNAs)),]
interactions <- interactions[which(interactions$effect %in% colnames(PE_matchedData$mRNAs)),]
nodes <- unique(union(interactions$cause, interactions$effect))

# TFs: Download the TF list from http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19
tfs <- read.csv(paste(rootDir, "/Data/Browse Transcription Factors hg19 - resource_browser.csv",
                      sep = ""))
i <- which(levels(tfs$Symbol) %in% nodes)
# extract TF exprssion data
tfData <- PE_matchedData$mRNAs[, levels(tfs$Symbol)[i]]

# Update PE data of mRNAs, only contain mRNAs in the network and not in the TF list
PE_matchedData$mRNAs <- PE_matchedData$mRNAs[,
                                             nodes[which(!(nodes %in% levels(tfs$Symbol)[i]))]]
mRNAsData_PE <-  PE_matchedData$mRNAs

# Combine data to a single matrix with columns are miRNA, mRNA, TFs
nomR <- ncol(PE_matchedData$mRNAs)
noTF <- ncol(tfData)
PE_data <- cbind(mRNAsData_PE, tfData)

# Free the memory
gc()

# Build the network
PE_network <-  buildNetworkForPersonalised(interactions, nomR, noTF, PE_data, 
                                           rootDir,usingpVal = TRUE, cutoff = 0.05)

# Save the network
write.csv(PE_network, paste(outDir, "/pVal_PE_network.csv", sep = ""), row.names = FALSE)

# Analyse network summary of the network
# PE_network <- read.csv(paste(outDir, "/PE_network.csv", sep = ""))

analyseNetworkForPersonalised(nomR, noTF, PE_network, PE_data, 
                              paste(outDir, "/pVal_PE_network_analysis.txt", sep = ""))

#================================================================
# (2) Identifying PE drivers
#================================================================

# Analyse controllability of the network
interactions <- read.csv(paste(outDir, "/pVal_PE_network.csv",
                               sep = ""))
# Write the edges of the network for analysing controllability
write.table(interactions, paste(outDir, "/Controllability/pVal_edges.dat", sep = ""), row.names = FALSE, col.names=FALSE, quote=FALSE)
# Run the controllability analysis
cmd <- paste(controlDir, "/parse.exe ", outDir, 
             "/Controllability/pVal_edges.dat", sep = "")
system(cmd)
cmd <- paste(controlDir, "/ControllabilityAnalysis.exe ", outDir, 
             "/Controllability/pVal_edges.dat", sep = "")
system(cmd)
# Analyse controllability of the network and output in a file
analyseControllability(paste(outDir, "/Controllability/pVal_edges.dat.output", sep = ""),
                       paste(outDir, "/pVal_analyseControllability.txt", sep = ""))

# Identify critical nodes in the network
# Read the result
nodetype <- read.table(paste(outDir, "/Controllability/pVal_edges.dat.nodetype", sep = ""))
colnames(nodetype) <- c("Name", "K", "Kin", "Kout", "TypeI", "TypeII")
nodetype$Name = as.character(nodetype$Name)
# Type 1 critical nodes of the network
Type1_critical_nodes <- nodetype[which(nodetype$TypeII == 0),]
Type1_critical_nodes$Name =  as.character(Type1_critical_nodes$Name)
Type2_critical_nodes <- nodetype[which(nodetype$TypeI == 0),]
Type2_critical_nodes$Name =  as.character(Type2_critical_nodes$Name)
# Save file
write.csv(Type1_critical_nodes, paste(outDir, "/pVal_Type1_critical_nodes.csv", sep = ""),
          row.names = FALSE)
write.csv(Type2_critical_nodes, paste(outDir, "/pVal_Type2_critical_nodes.csv", sep = ""),
          row.names = FALSE)

# Combine Control data, rank candidate PE drivers
index = c(which(pd$subtype == type),which(pd$subtype == "Control" & pd$GA >= 34 & pd$GA != 'NA'))
design <- model.matrix(~0+factor(pd$subtype[index]))
colnames(design)=levels(factor(pd$subtype[index]))
rownames(design)=colnames(ddata[,index])
design

index = c(which(pd$subtype == type),which(pd$phenotype == "normal" & pd$GA >= 34))

temp = pd$subtype[index]
temp[temp=="NA"] <- "Control"
design <- model.matrix(~0+factor(temp))
colnames(design)=levels(factor(temp))
rownames(design)=colnames(ddata[,index])
design

contrast.matrix<-makeContrasts(paste0(unique(temp),collapse = "-"),levels = design)
contrast.matrix 

contrast.matrix<-makeContrasts(paste0(unique(pd$subtype[index]),collapse = "-"),levels = design)
contrast.matrix 

##step1
fit <- lmFit(ddata[,index],design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
limma_notrend_results_LoPE = na.omit(tempOutput) 
head(limma_notrend_results_LoPE)
write.csv(limma_notrend_results_LoPE, paste0(outDir, "/limma_notrend_results_LoPE.csv"),quote = F)

limma_notrend_results_LoPE$gene = rownames(limma_notrend_results_LoPE)
limma_notrend_results_LoPE = limma_notrend_results_LoPE[which(limma_notrend_results_LoPE$adj.P.Val < 0.05),c(7,1)]
limma_notrend_results_LoPE$logFC = abs(limma_notrend_results_LoPE$logFC)

Type1_candidate_LoPE_drivers <- merge(Type1_critical_nodes, limma_notrend_results,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")


Type2_candidate_LoPE_drivers <- merge(Type2_critical_nodes, limma_notrend_results,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")

# By the absolute log-fold-change
Type1_coding_candidate_LoPE_drivers <- 
  Type1_candidate_LoPE_drivers[order(Type1_candidate_LoPE_drivers$logFC, decreasing = TRUE),]

Type2_coding_candidate_LoPE_drivers <- 
  Type2_candidate_LoPE_drivers[order(Type2_candidate_LoPE_drivers$logFC, decreasing = TRUE),]


# Compare EoPE and LoPE biomarkers
length(which(! Type1_coding_candidate_LoPE_drivers$Name %in% Type1_candidate_EoPE_drivers$Name)) #LoPE
length(which(! Type1_candidate_EoPE_drivers$Name %in% Type1_coding_candidate_LoPE_drivers$Name)) #EoPE
length(which(Type1_coding_candidate_LoPE_drivers$Name %in% Type1_candidate_EoPE_drivers$Name)) #both
length(which(Type2_coding_candidate_LoPE_drivers$Name %in% Type2_candidate_EoPE_drivers$Name)) #bot

index = which(! Type1_coding_candidate_LoPE_drivers$Name %in% Type1_candidate_EoPE_drivers$Name)
type1_LoPE = Type1_coding_candidate_LoPE_drivers[index,] #265
index = which(alias2SymbolTable(type1_LoPE$Name) %in% alias2Symbol(gold_standard$dbPEC))
type1_LoPE$dbPEC = "No"
type1_LoPE$dbPEC[index] = "Yes"

index = which(Type1_candidate_EoPE_drivers$Name %in% Type1_coding_candidate_LoPE_drivers$Name)
type1_EoPE = Type1_candidate_EoPE_drivers[-index,] #298
index = which(alias2SymbolTable(type1_EoPE$Name) %in% alias2Symbol(gold_standard$dbPEC))
type1_EoPE$dbPEC = "No"
type1_EoPE$dbPEC[index] = "Yes"

index = which(Type1_coding_candidate_LoPE_drivers$Name %in% Type1_candidate_EoPE_drivers$Name)
type1_both = Type1_coding_candidate_LoPE_drivers[index,] #418
index = which(alias2SymbolTable(type1_both$Name) %in% alias2Symbol(gold_standard$dbPEC))
type1_both$dbPEC = "No"
type1_both$dbPEC[index] = "Yes"

index = which(Type2_coding_candidate_LoPE_drivers$Name %in% Type2_candidate_EoPE_drivers$Name)
Type2_LoPE = Type2_coding_candidate_LoPE_drivers[-index,] #469
type2_both = Type2_coding_candidate_LoPE_drivers[index,] #563

index = which(alias2SymbolTable(Type2_LoPE$Name) %in% alias2Symbol(gold_standard$dbPEC))
Type2_LoPE$dbPEC = "No"
Type2_LoPE$dbPEC[index] = "Yes"
index = which(alias2SymbolTable(type2_both$Name) %in% alias2Symbol(gold_standard$dbPEC))
type2_both$dbPEC = "No"
type2_both$dbPEC[index] = "Yes"

index = which(Type2_candidate_EoPE_drivers$Name %in% Type2_coding_candidate_LoPE_drivers$Name)
Type2_EoPE = Type2_candidate_EoPE_drivers[-index,] #494
index = which(alias2SymbolTable(Type2_EoPE$Name) %in% alias2Symbol(gold_standard$dbPEC))
Type2_EoPE$dbPEC = "No"
Type2_EoPE$dbPEC[index] = "Yes"


#write to a xlsx file
library("xlsx")
write.xlsx(type1_EoPE, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type2_EoPE", append = TRUE)
write.xlsx(type1_LoPE, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type2_LoPE", append = TRUE)
write.xlsx(type1_both, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type2_both", append = TRUE)
write.xlsx(Type2_EoPE, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type1_EoPE", append = TRUE)
write.xlsx(Type2_LoPE, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type1_LoPE", append = TRUE)
write.xlsx(type2_both, file = "subtypes2.xlsx", row.names = F,
           sheetName = "Type1_both", append = TRUE)
