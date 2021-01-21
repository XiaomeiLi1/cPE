## Data dived PPI but filted by  the directed PPI

# This is the script of the proposed method and its application in detecting preeclampsia biomarkers. 
# The method is applied to the preeclampsia dataset to discover preeclampsia biomarkers.
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
# PE biomarkers
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
outDir <- "C:/unisa/RstudioProjects/PEDriver/result/PE" # Output folder
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
PE_matchedData$PEType = "preeclampsia"
index = which(pd$phenotype == "preeclampsia")
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
PE_data <- cbind(miRNAsData_PE, mRNAsData_PE, tfData)

# Free the memory
gc()

# Build the network
PE_network <-  buildNetworkForPersonalised(interactions, nomR, noTF, PE_data, rootDir,usingpVal = FALSE, cutoff = 0.05)

# Save the network
write.csv(PE_network, paste(outDir, "/pVal_PE_network.csv", sep = ""), row.names = FALSE)

# Analyse network summary of the network
# PE_network <- read.csv(paste(outDir, "/PE_network.csv", sep = ""))

analyseNetworkForPersonalised(nomR, noTF, PE_network, PE_data, 
                              paste(outDir, "/pVal_PE_network_analysis.txt", sep = ""))

#================================================================
# (2) Identifying PE drivers
#================================================================
# Load the mutation data
PESNPdb = read.csv("C:/unisa/RstudioProjects/PEDriver/Data/RAW/PESNPdb.csv", header = T, stringsAsFactors = F)
temp = PESNPdb[,c(3,2)]
temp <- unique(temp)
temp = table(temp[,1])
proteinAffectingMut = data.frame(row.names = names(temp))
proteinAffectingMut$symbol = names(temp)
proteinAffectingMut$weight  = as.vector(temp)

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


# plot of degrees
avg_din = rep(0,3)
avg_din[1] = mean(nodetype[which(nodetype$TypeII == 0),"Kin"])
avg_din[3] = mean(nodetype[which(nodetype$TypeII == 1),"Kin"])
avg_din[2] = mean(nodetype[which(nodetype$TypeII == 2),"Kin"])
names(avg_din) = c("Critical", "Ordinary", "Redundant")

# Type 1 in degree
df=data.frame(row.names = names(avg_din))
df$x = names(avg_din)
df$x = factor(df$x, level=names(rev(avg_din)))
df$y = avg_din

p1=ggplot(data=df, aes(x=x, y=y,fill = x)) +
  xlab("") + 
  ylab("Average in-degree") + 
  #labs(tag = "A")+
  scale_fill_manual(values=c("gray30", "Blue3", "Red3"))+
  geom_bar(stat="identity",width=0.5,color="black")+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,6),
                     breaks  = seq(0,6, by = 1)) +
  coord_flip()

avg_dout = rep(0,3)
avg_dout[1] = mean(nodetype[which(nodetype$TypeII == 0),"Kout"])
avg_dout[3] = mean(nodetype[which(nodetype$TypeII == 1),"Kout"])
avg_dout[2] = mean(nodetype[which(nodetype$TypeII == 2),"Kout"])
names(avg_dout) = c("Critical", "Ordinary", "Redundant")

#Out degress
df=data.frame(row.names = names(avg_dout))
df$x = names(avg_dout)
df$x = factor(df$x, level=names(rev(avg_dout)))
df$y = avg_dout

p2=ggplot(data=df, aes(x=x, y=y,fill = x)) +
  xlab("") + 
  ylab("Average out-degree") + 
  #labs(tag = "B")+
  scale_fill_manual(values=c("gray30", "Blue3", "Red3"))+
  geom_bar(stat="identity",width=0.5,color="black")+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,4),
                     breaks  = seq(0,4, by = 1)) +
  coord_flip()

# Type 2 in degree
avg_dout = rep(0,3)
avg_din[1] = mean(nodetype[which(nodetype$TypeI == 0),"Kin"])
avg_din[3] = mean(nodetype[which(nodetype$TypeI == 1),"Kin"])
avg_din[2] = mean(nodetype[which(nodetype$TypeI == 2),"Kin"])
names(avg_din) = c("Critical", "Ordinary", "Redundant")

df=data.frame(row.names = names(avg_din))
df$x = names(avg_din)
df$x = factor(df$x, level=names(rev(avg_din)))
df$y = avg_din

p11=ggplot(data=df, aes(x=x, y=y,fill = x)) +
  xlab("") + 
  ylab("Average in-degree") + 
  #labs(tag = "A")+
  scale_fill_manual(values=c("gray30", "Blue3", "Red3"))+
  geom_bar(stat="identity",width=0.5,color="black")+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,8),
                     breaks  = seq(0,8, by = 2)) +
  coord_flip()

avg_dout = rep(0,3)
avg_dout[1] = mean(nodetype[which(nodetype$TypeI == 0),"Kout"])
avg_dout[3] = mean(nodetype[which(nodetype$TypeI == 1),"Kout"])
avg_dout[2] = mean(nodetype[which(nodetype$TypeI == 2),"Kout"])
names(avg_dout) = c("Critical", "Ordinary", "Redundant")

df=data.frame(row.names = names(avg_dout))
df$x = names(avg_dout)
df$x = factor(df$x, level=names(rev(avg_dout)))
df$y = avg_dout

p22=ggplot(data=df, aes(x=x, y=y,fill = x)) +
  xlab("") + 
  ylab("Average out-degree") + 
  #labs(tag = "B")+
  scale_fill_manual(values=c("gray30", "Blue3", "Red3"))+
  geom_bar(stat="identity",width=0.5,color="black")+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,8),
                     breaks  = seq(0,8, by = 2)) +
  coord_flip()


#Accumulative in-degree distribution
library(igraph)

Binr <- read.csv(paste(outDir, "/pVal_PE_network.csv", sep = ""), header = T, stringsAsFactors = F)
g <- graph.data.frame(Binr,directed = TRUE)
#din <- degree(g, v = Type1_critical_nodes$Name, mode="in")
ddinc <- degree.distribution(g, v = Type1_critical_nodes$Name, mode="in", cumulative=TRUE)
ddino <- degree.distribution(g, v = nodetype[which(nodetype$TypeI == 2),"Name"], mode="in", cumulative=TRUE)
ddinr <- degree.distribution(g, v = nodetype[which(nodetype$TypeI == 1),"Name"], mode="in", cumulative=TRUE)

ddinc <- degree.distribution(g, v = Type2_critical_nodes$Name, mode="in", cumulative=TRUE)
ddino <- degree.distribution(g, v = nodetype[which(nodetype$TypeII == 2),"Name"], mode="in", cumulative=TRUE)
ddinr <- degree.distribution(g, v = nodetype[which(nodetype$TypeII == 1),"Name"], mode="in", cumulative=TRUE)

df=data.frame(row.names = 1:(length(ddinc)+length(ddino)+length(ddinr)))
df$x = c(1:length(ddinc),1:length(ddino),1:length(ddinr))
df$y = c(ddinc,ddino,ddinr)
df$group = c(rep("Critical",length(ddinc)),rep("Ordinary",length(ddino)),rep("Redundant",length(ddinr)))

p3 <- ggplot(df,aes(x=x,y=y,color=group))+
  scale_color_manual(values = c("red3", "blue3", "gray30")) +
  xlab("In-degree") + 
  ylab("Probability")+
  labs(tag = "C")+
  geom_point(size = 3,shape=16) +
  theme(#plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.title = element_blank(),
    legend.position="bottom",
    #legend.box = "horizontal",
    #legend.box.background = element_rect(colour = "black",linetype="solid"),
    legend.text = element_text(face = "bold", size = 10))+
  coord_trans(x="log10", y="log10")


ddoutc <- degree.distribution(g, v = Type1_critical_nodes$Name, mode="out", cumulative=TRUE)
ddouto <- degree.distribution(g, v = nodetype[which(nodetype$TypeI == 2),"Name"], mode="out", cumulative=TRUE)
ddoutr <- degree.distribution(g, v = nodetype[which(nodetype$TypeI == 1),"Name"], mode="out", cumulative=TRUE)

ddoutc <- degree.distribution(g, v = Type2_critical_nodes$Name, mode="out", cumulative=TRUE)
ddouto <- degree.distribution(g, v = nodetype[which(nodetype$TypeII == 2),"Name"], mode="out", cumulative=TRUE)
ddoutr <- degree.distribution(g, v = nodetype[which(nodetype$TypeII == 1),"Name"], mode="out", cumulative=TRUE)

df=data.frame(row.names = 1:(length(ddoutc)+length(ddouto)+length(ddoutr)))
df$x = c(1:length(ddoutc),1:length(ddouto),1:length(ddoutr))
df$y = c(ddoutc,ddouto,ddoutr)
df$group = c(rep("Critical",length(ddoutc)),rep("Ordinary",length(ddouto)),rep("Redundant",length(ddoutr)))

p4 <- ggplot(df,aes(x=x,y=y,color=group))+
  scale_color_manual(values = c("red3", "blue3", "gray30")) +
  xlab("Out-degree") + 
  ylab("Probability")+
  labs(tag = "D")+
  geom_point(size = 3,shape=16) +
  theme(#plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.title = element_blank(),
    legend.position="bottom",
    #legend.box = "horizontal",
    #legend.box.background = element_rect(colour = "black",linetype="solid"),
    legend.text = element_text(face = "bold", size = 10))+
  coord_trans(x="log10", y="log10")

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggrepel)

pdf("C:/unisa/RstudioProjects/PEDriver/Figure1-1.pdf",width=7,height=6,onefile = FALSE)

grid.arrange(p1, p2, p3, p4,
             ncol = 2, heights = 1:2)


dev.off()

pdf("C:/unisa/RstudioProjects/PEDriver/Figure2-1.pdf",width=7,height=6,onefile = FALSE)
grid.arrange(p11, p22, p3, p4,
             ncol = 2, heights = 1:2)
dev.off()

# DE analysis
# Combine Control data, rank candidate PE drivers
index = c(which(pd$phenotype == "preeclampsia"),which(pd$phenotype == "normal" & pd$GA >= 34))
design <- model.matrix(~0+factor(pd$phenotype[index]))
colnames(design)=levels(factor(pd$phenotype[index]))

ddata = exprs(GSE75010)

rownames(design)=colnames(ddata[,index])
design

contrast.matrix<-makeContrasts(paste0(unique(pd$phenotype[index]),collapse = "-"),levels = design)
contrast.matrix 

##step1
fit <- lmFit(ddata[,index],design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
limma_notrend_results = na.omit(tempOutput) 
head(limma_notrend_results)
write.csv(limma_notrend_results, paste0(outDir, "/limma_notrend_results_2.csv"),quote = F)

limma_notrend_results = read.csv("C:/unisa/RstudioProjects/PEDriver/result/PE/limma_notrend_results_2.csv",header = T, stringsAsFactors = F)
rownames(limma_notrend_results) = limma_notrend_results$X
colnames(limma_notrend_results)[1]="gene"
DE_limma = limma_notrend_results[nodes,]
DE_critical = DE_limma[which(limma_notrend_results$adj.P.Val < 0.05),1:2]
limma_notrend_results = limma_notrend_results[which(limma_notrend_results$adj.P.Val < 0.05),1:2]
limma_notrend_results$logFC = abs(limma_notrend_results$logFC)

Type1_candidate_PE_drivers <- merge(Type1_critical_nodes, limma_notrend_results,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")

Type2_candidate_PE_drivers <- merge(Type2_critical_nodes, limma_notrend_results,
                                    all.x = TRUE, by.x = "Name", by.y = "gene")


# By the absolute log-fold-change
Type1_coding_candidate_PE_drivers <- 
  Type1_candidate_PE_drivers[order(Type1_candidate_PE_drivers$logFC, decreasing = TRUE),]

Type2_coding_candidate_PE_drivers <- 
  Type2_candidate_PE_drivers[order(Type2_candidate_PE_drivers$logFC, decreasing = TRUE),]

# Comparison methods
# DE
DE_ranking = limma_notrend_results[order(limma_notrend_results$logFC, decreasing = TRUE),]
critical_nodes = rbind(Type1_candidate_PE_drivers,Type2_candidate_PE_drivers)
critical_nodes = critical_nodes[order(critical_nodes$logFC, decreasing = TRUE),]

# characterize top 20 genes of critical nodes
Top50 = data.frame(row.names = 1:50)
Top50$gene = critical_nodes$Name[1:50]
Top50$Type = rep(1,50)
index = which(critical_nodes$TypeII[1:50] == 0)
Top50$Type[index] = 2
index = which(alias2Symbol(Top50$gene)%in%alias2Symbol(gold_standard$dbPEC))
Top50$dbPEC = rep("No", 50)
Top50$dbPEC[index] = "Yes"
write.csv(Top50, file = "Top50.csv")


# Hubs
library(igraph)
g1 <- graph.data.frame(interactions, directed=TRUE)
hub=hub_score(g1, scale = TRUE, weights = NULL, options = arpack_defaults)
hub_rank = hub$vector[order(hub$vector,decreasing=TRUE) [which(hub$vector > 0)]]


#================================================================
# Validation
#================================================================
# Load the gold standard, CGC


df=data.frame(row.names = 1:10)
df$x = c("Top 50","Top 100","Top 150","Top 200","All","Top 50","Top 100","Top 150","Top 200","All")
df$x = factor(df$x, levels = c("Top 50","Top 100","Top 150","Top 200","All"))
df$y = c(1:10)
df$group = c(rep("Type-2",5),rep("Type-1",5))

df=data.frame(row.names = 1:6)
df$x = c("Top 50","Top 200","All","Top 50","Top 200","All")
df$x = factor(df$x, levels = c("Top 50","Top 200","All"))
df$y = c(1:6)
df$group = c(rep("Type-2",3),rep("Type-1",3))

df$y[1] = length(which(alias2SymbolTable(Type1_coding_candidate_PE_drivers$Name[1:50])%in%alias2Symbol(gold_standard$dbPEC)))
df$y[2] = length(which(alias2SymbolTable(Type1_coding_candidate_PE_drivers$Name[1:200])%in%alias2Symbol(gold_standard$dbPEC)))
df$y[3] = length(which(alias2SymbolTable(Type1_coding_candidate_PE_drivers$Name)%in%alias2Symbol(gold_standard$dbPEC)))
df$y[4] = length(which(alias2SymbolTable(Type2_coding_candidate_PE_drivers$Name[1:50])%in%alias2Symbol(gold_standard$dbPEC)))
df$y[5] = length(which(alias2SymbolTable(Type2_coding_candidate_PE_drivers$Name[1:200])%in%alias2Symbol(gold_standard$dbPEC)))
df$y[6] = length(which(alias2SymbolTable(Type2_coding_candidate_PE_drivers$Name)%in%alias2Symbol(gold_standard$dbPEC)))


p5 = ggplot(data=df, aes(x=x, y=y,fill = group)) +
  xlab("Top N genes") + 
  ylab("Number of validated genes") + 
  labs(fill = "Node type")+
  scale_fill_manual(values=c("steelblue", "tomato"))+
  geom_bar(stat="identity",width=0.5,color="black", position=position_dodge())+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position =  c(0.2,0.7),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,255),
                     breaks  = seq(0,250, by = 50))+
  geom_text(aes(label=y), position=position_dodge(width=0.5), vjust=-0.25)
p5

##Comparison methods
dff = data.frame(row.names = c("DE", "Hub", "cPE"))
dff$x = c("DE", "Hub", "cPE")
dff$x = factor(dff$x, levels = rownames(dff))
dff$y = c(831.0/4096,998.0/4370, 429.0/1722)

p6 = ggplot(data=dff, aes(x=x, y=y)) +
  xlab("Method") + 
  ylab("Percentage of validated genes") + 
  geom_bar(stat="identity",width=0.5,color="black", fill = "grey30", position=position_dodge())+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #panel.background = element_blank(),
    #panel.grid.major = element_blank(),
    legend.position =  "none",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,0.3),
                     breaks  = seq(0,0.3, by = 0.1),
                     labels = scales::percent_format(accuracy = 1))+
  geom_text(aes(label=scales::percent(y)), position=position_dodge(width=0.9), vjust=-0.25)
p6

pdf("C:/unisa/RstudioProjects/PEDriver/Validation.pdf",width=7,height=4,onefile = FALSE)
ggarrange(p5, p6, widths = c(4,3))
dev.off()

# This function allows you to get drivers which are specific for cPE
setdiff(critical_nodes$Name, c(DE_ranking$gene, names(hub_rank)))

# P value of cPE
pValue(429,1138, 5489-1138, 1722) #1.962691e-07

# Overlap between different methods
library(VennDiagram)

venn.plot.top50 <- venn.diagram(
  x = list(DE_ranking_2$gene[1:50], names(hub_rank)[1:50],critical_nodes_2$Name[1:50]),
  category.names = c("DE" , "Hub" , "cPE"),
  filename = NULL,
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3))
)

venn.plot.top200 <- venn.diagram(
  x = list(DE_ranking_2$gene[1:200], names(hub_rank)[1:200],critical_nodes_2$Name[1:200]),
  category.names = c("DE" , "Hub" , "cPE"),
  filename = NULL,
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3))
)

pdf("Venn2.pdf",width=7,height=3.5,onefile = FALSE)
grid.arrange(gTree(children=venn.plot.top50), gTree(children=venn.plot.top200), ncol=2)
dev.off()

### Enrichment analysis
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

index = which(alias2SymbolTable(critical_nodes$Name)%in%alias2Symbol(gold_standard$dbPEC))
novel = critical_nodes$Name[-index]

novel_GO = bitr(novel, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
novel_ALL <- enrichGO(gene = novel_GO$ENTREZID, 
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, 
                    readable = TRUE) 

novel_ALL_table = novel_ALL@result

dbPEC_GO = bitr(alias2SymbolTable(gold_standard$dbPEC), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
dbPEC_ALL <- enrichGO(gene = dbPEC_GO$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, 
                      readable = TRUE) 

dbPEC_ALL_table = dbPEC_ALL@result

index = which(novel_ALL_table$ID %in% dbPEC_ALL_table$ID)

GO_stat = data.frame(row.names = 1:nrow(novel_ALL_table))
GO_stat = novel_ALL_table[order(novel_ALL_table$p.adjust, decreasing = FALSE),c(1,2,3,7,9)]
GO_stat = novel_ALL_table[order(novel_ALL_table$p.adjust, decreasing = FALSE),]
GO_stat = GO_stat[-index,]
GO_stat$dbPEC = rep("No", nrow(novel_ALL_table))
GO_stat$dbPEC[index] = "Yes"
write.csv(GO_stat, file = "GO.csv", row.names = F)

ego_ALL = novel_ALL
ego_ALL@result = GO_stat
pdf("C:/unisa/RstudioProjects/PEDriver/Figure6.pdf",width=10,height=6,onefile = FALSE)
barplot(ego_ALL, showCategory=20,title="Enrichment GO")
dev.off()

pdf("C:/unisa/RstudioProjects/PEDriver/GOALL.pdf",width=10,height=6,onefile = FALSE)
barplot(novel_ALL, showCategory=20,title="Enrichment GO")
dev.off()