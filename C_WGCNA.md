# WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA)

## Re-organising and filtering data
```{R}
## load libraries
rm(list=ls())
library(WGCNA)
library(cluster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(flashClust)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## filtering genes
refGenes <- scan("keep_genes_protein_coding.txt", what = "", sep = "\n") # 21578
genesMito <- scan("mtGenes.txt", what = "", sep = "\n") # 13
keepGenes <- setdiff(refGenes, genesMito) # 21565

## load counts, rpkm, cpm table
load("VS_COUNTS_RPKM_CPM_META.RData")

## load updated meta data 
newmeta <- read.table("VS_META_TABLE2.txt", sep = "\t", header = TRUE)
row.names(newmeta) <- newmeta$SAMPLE_1
newmeta <- newmeta[,c("GENOTYPE", "AGE", "SEX", "RNA_EXT_BATCH", "HARVEST_BATCH", "RIN", "SAMPLE")]

## re-rganize counts, rpkm and cpm tables
newmeta$GENOAGESAMPLE <- paste(newmeta$GENOTYPE, newmeta$AGE, newmeta$SAMPLE, sep = "_")
meta2 <- as.data.frame(newmeta[,c("SAMPLE", "GENOAGESAMPLE", "GENOTYPE", "AGE")])
counts2 <- merge(as.data.frame(t(counts)), meta2, by = "row.names")
rpkm2 <- merge(as.data.frame(t(rpkm)), meta2, by = "row.names")
cpm2 <- merge(as.data.frame(t(cpm)), meta2, by = "row.names")

row.names(newmeta) <- newmeta$GENOAGESAMPLE
row.names(counts2) <- counts2$GENOAGESAMPLE
row.names(rpkm2) <- rpkm2$GENOAGESAMPLE
row.names(cpm2) <- cpm2$GENOAGESAMPLE

counts2$Row.names <- NULL
counts2$SAMPLE <- NULL
counts2$GENOAGESAMPLE <- NULL
counts2$GENOTYPE <- NULL
counts2$AGE <- NULL

rpkm2$Row.names <- NULL
rpkm2$SAMPLE <- NULL
rpkm2$GENOAGESAMPLE <- NULL
rpkm2$GENOTYPE <- NULL
rpkm2$AGE <- NULL

cpm2$Row.names <- NULL
cpm2$SAMPLE <- NULL
cpm2$GENOAGESAMPLE <- NULL
cpm2$GENOTYPE <- NULL
cpm2$AGE <- NULL

DESIGN <- newmeta[order(row.names(newmeta)),]
COUNTS <- as.data.frame(t(counts2[row.names(counts2) %in% row.names(DESIGN),colnames(counts2) %in% keepGenes]))
COUNTS <- COUNTS[,order(colnames(COUNTS))]
RPKM <- as.data.frame(t(rpkm2[row.names(rpkm2) %in% row.names(DESIGN),colnames(rpkm2) %in% keepGenes]))
RPKM <- RPKM[,order(colnames(RPKM))]
CPM <- as.data.frame(t(cpm2[row.names(cpm2) %in% row.names(DESIGN),colnames(cpm2) %in% keepGenes]))
CPM <- CPM[,order(colnames(CPM))]

## making sure the row names in colData matches to column names in counts_data
all(colnames(COUNTS) %in% rownames(DESIGN))
all(colnames(RPKM) %in% rownames(DESIGN))
all(colnames(CPM) %in% rownames(DESIGN))

## making sure they in the same order
all(colnames(COUNTS) == rownames(DESIGN))
all(colnames(RPKM) == rownames(DESIGN))
all(colnames(CPM) == rownames(DESIGN))

## remove outliers
outliers6M <- c("WT_06m_SCA6_7_4", "WT_06m_SCA6_8_5", "HET_06m_SCA6_7_8", "HET_06m_SCA6_7_9")
outliers12M <- c("WT_12m_SCA6_2_7", "HET_12m_SCA6_3_1")
outliers <- c(outliers6M, outliers12M)

myDesign <- DESIGN[!row.names(DESIGN) %in% outliers,]
myCOUNTS <- COUNTS[,colnames(COUNTS) %in% row.names(myDesign)]
myCPM <- CPM[,colnames(CPM) %in% row.names(myDesign)]
myRPKM <- RPKM[,colnames(RPKM) %in% row.names(myDesign)]

## filtering genes (CPM >= 1 across GENOTYPEA_AGE)
filter <- apply(myCPM, 1, function(x) { all(x[1:5] >= 1) | all(x[6:8] >= 1) | all(x[9:12] >= 1) | all(x[13:17] >= 1) | all(x[18:22] >= 1) | all(x[23:25] >= 1) | all(x[26:29] >= 1) | all(x[30:34] >= 1)})
CPMfilt <- myCPM[filter,]
COUNTSfilt <- myCOUNTS[filter,]
RPKMfilt <- myRPKM[filter,]

## update meta or design
myDesign$BATCH <- paste(myDesign$RNA_EXT_BATCH, myDesign$HARVEST_BATCH, sep = "_")

tab <- log2(RPKMfilt + 1)
tab[is.na(tab)] <- 0

datExpr <- as.data.frame(t(tab))
names(datExpr) <- rownames(tab)
rownames(datExpr) <- names(tab)

datExpr0 <- datExpr
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

## filter out the missing entries from the data
if (!gsg$allOK)
 {
 if (sum(!gsg$goodGenes) > 0)
  printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
 if (sum(!gsg$goodSamples) > 0)
  printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
 datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
 }

## perform hirerchical clustering on samples (not genes)
sampleTree <- hclust(dist(datExpr), method="average")

pdf("VS_sampleTree1.pdf", width=8, height=6)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## re-arranging traits data
traitsData <- myDesign[,c("GENOTYPE", "AGE")]
datTraits <- as.data.frame(cbind(as.numeric(factor(traitsData$GENOTYPE)), as.numeric(factor(traitsData$AGE))))
colnames(datTraits) <- c("GENOTYPE", "AGE")
row.names(datTraits) <- rownames(traitsData)

## perform hirerchical clustering on samples (not genes) along with heatmap for traits
sampleTree2 <- hclust(dist(datExpr), method="average")
traitColors <- labels2colors(traitsData[rownames(datExpr),])

pdf("VS_sampleTree2.pdf", width=8, height=6)
par(cex=0.6)
par(mar=c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), main="Sample dendrogram and trait heatmap")
#abline(h=15, col="red")
dev.off()

## save the processed data to a R file
save(datExpr, datTraits, file="VS_WGCNA_DATA.RData")
```

&nbsp;

## Network Construction and Modules Detection
```{R}
## soft thresholding
powers <- c(seq(2, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = 20000, networkType = "signed") 

#  pickSoftThreshold: calculating connectivity for given powers...
#    ..working on genes 1 through 13374 of 13374
#    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      2   0.1300  2.850          0.946  3870.0   3850.00 4820.0
# 2      4   0.0927 -0.988          0.968  1530.0   1490.00 2490.0
# 3      6   0.3070 -1.390          0.975   722.0    691.00 1550.0
# 4      8   0.4920 -1.570          0.977   386.0    358.00 1040.0
# 5     10   0.6190 -1.690          0.979   225.0    200.00  747.0
# 6     12   0.7150 -1.790          0.988   140.0    118.00  555.0
# 7     14   0.7750 -1.850          0.992    91.7     72.80  425.0
# 8     16   0.8140 -1.890          0.995    62.5     46.60  333.0
# 9     18   0.8390 -1.900          0.994    44.0     30.60  265.0
# 10    20   0.8550 -1.920          0.994    31.9     20.60  217.0
# 11    22   0.8580 -1.960          0.994    23.6     14.20  180.0
# 12    24   0.8660 -1.980          0.995    17.9      9.92  152.0
# 13    26   0.8750 -1.970          0.996    13.7      7.06  129.0
# 14    28   0.8800 -1.970          0.996    10.7      5.09  110.0
# 15    30   0.8830 -1.970          0.993     8.5      3.73   95.4


pdf("VS_SoftThresholdingPower_signed.pdf", width=8, height=4)
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue"); abline(h=0.9,col="green"); 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

## network construction
PWR <- 16
net <- blockwiseModules(
       # Input data
       datExpr, 
       # Data checking options
       checkMissingData = FALSE,
       # Options for splitting data into blocks
       #blocks = NULL,
       maxBlockSize = 10000,
       #blockSizePenaltyPower = 5,
       #nPreclusteringCenters = as.integer(min(ncol(datExpr)/20, 100*ncol(datExpr)/maxBlockSize)),
       #randomSeed = 12345,
       # load TOM from previously saved file?
       #loadTOM = FALSE,
       # Network construction arguments: correlation options
       corType = "pearson",
       #maxPOutliers = 1, 
       #quickCor = 0,
       #pearsonFallback = "individual",
       #cosineCorrelation = FALSE,
       # Adjacency function options
       power = PWR,
       networkType = "signed",
       #replaceMissingAdjacencies = FALSE,
       # Topological overlap options
       TOMType = "signed",
       TOMDenom = "mean",
       # Saving or returning TOM
       #getTOMs = NULL,
       saveTOMs = TRUE, 
       saveTOMFileBase = "VS_Net_blockwiseTOM",
       # Basic tree cut options
       # deepSplit = 4,
       # detectCutHeight = 0.995, 
       minModuleSize = 50, # min(20, ncol(datExpr)/2 ),
       # Advanced tree cut options
       #maxCoreScatter = NULL, minGap = NULL,
       #maxAbsCoreScatter = NULL, minAbsGap = NULL,
       #minSplitHeight = NULL, minAbsSplitHeight = NULL,
       #useBranchEigennodeDissim = FALSE,
       #minBranchEigennodeDissim = mergeCutHeight,
       #stabilityLabels = NULL,
       #minStabilityDissim = NULL,
       pamStage = TRUE, 
       pamRespectsDendro = TRUE,
       # Gene reassignment, module trimming, and module "significance" criteria
       # reassignThreshold = 1e-2,
       #minCoreKME = 0.5, 
       #minCoreKMESize = minModuleSize/3,
       #minKMEtoStay = 0.3,
       # Module merging options
       # mergeCutHeight = 0.05, 
       #impute = TRUE, 
       #trapErrors = FALSE, 
       # Output options
       numericLabels = TRUE,
       # Options controlling behaviour
       nThreads = 48,
       verbose = 0#, indent = 0
       )

## detected modules and size of each module (number of genes in each module)
## Module 0 is the module containing genes which do not fall into any of the module
table(net$colors)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1980  969  965  704  531  471  459  448  404  367  362  321  313  306  275  263 
#   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
#  253  218  216  209  203  197  196  171  170  169  166  154  153  151  147  140 
#   32   33   34   35   36   37   38   39   40   41   42   43   44 
#  137  133  131  123  113  112  101   89   84   83   82   68   67

moduleLabelsAutomatic <- net$colors
moduleColorsAutomatic <- labels2colors(moduleLabelsAutomatic)
MEsAutomatic <- net$MEs
unique(moduleColorsAutomatic)
#  [1] "tan"             "pink"            "yellow"          "yellowgreen"    
#  [5] "mediumpurple3"   "darkred"         "blue"            "darkgreen"      
#  [9] "lightyellow"     "turquoise"       "black"           "lightgreen"     
# [13] "paleturquoise"   "grey60"          "lightcyan"       "darkolivegreen" 
# [17] "plum1"           "grey"            "purple"          "royalblue"      
# [21] "cyan"            "skyblue"         "lightsteelblue1" "brown"          
# [25] "green"           "saddlebrown"     "darkturquoise"   "floralwhite"    
# [29] "red"             "steelblue"       "violet"          "salmon"         
# [33] "greenyellow"     "magenta"         "orange"          "darkgrey"       
# [37] "skyblue3"        "orangered4"      "lightcyan1"      "darkmagenta"    
# [41] "midnightblue"    "sienna3"         "darkorange"      "white"          
# [45] "ivory"

## save network
save(net, file="VS_Net_Signed_Pwr16_MinMidSize50.RData")
write.table(moduleColorsAutomatic, "VS_Net_Signed_Pwr16_MinMidSize50_Colors.txt", sep="\t", quote=F)


## plotting a dendrogram showing hirerchical clustering with each module indicated by its assigned color
mergedColors <- labels2colors(net$colors)
pdf("VS_Net_ModuleDendrogram.pdf", width=8, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

## KME
KMEs <- signedKME(datExpr, net$MEs, corFnc = "cor", corOptions = "use = 'p'")
kme <- data.frame(names(datExpr), moduleColorsAutomatic, KMEs)
kme$names.datExpr. <- NULL
write.table(kme, "VS_Net_Signed_Pwr16_MinMidSize50_KMEs.txt", sep="\t", quote=F)

# reclusterd without the grey
restGenes <- (moduleColorsAutomatic != "grey")
diss1 <- 1 - TOMsimilarityFromExpr(datExpr[,restGenes], power = PWR, corType = "pearson", networkType="signed", TOMType = "signed", TOMDenom = "mean", nThreads = 24, verbose = 5, indent = 0)
colnames(diss1) <- rownames(diss1) <- names(datExpr)[restGenes]
hier1 <- flashClust(as.dist(diss1), method="average" )
plotColors <- cbind(moduleColorsAutomatic)
plotColors <- as.data.frame(plotColors)
nogrey <- plotColors[plotColors$module != "grey",]

pdf("VS_Net_ModuleDendrogram_nogrey.pdf")
par(mfrow=c(2,1))
plotDendroAndColors(hier1, nogrey, dendroLabels = FALSE )
dev.off()

# Make Ggplot2 input for WGCNA eigengene values
# nGenes <- ncol(datExpr) 
# nSamples <- nrow(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColorsAutomatic, softPower = PWR, impute = TRUE)$eigengenes
MEs0$Rows <- colnames(tab[,row.names(datExpr)])
MEs0$Class <- c(rep("HET_03m", 5), rep("HET_06m", 3), rep("HET_12m", 4), rep("HET_19m", 5), rep("WT_03m", 5), rep("WT_06m", 3), rep("WT_12m", 4), rep("HET_19m", 5))
write.table(MEs0, "VS_Net_Matrix_Module_Correlation.txt",sep="\t",quote=F)

## define the association between genes and detected modules.
Adj <- adjacency(datExpr, power = PWR, type = "signed",corFnc = "cor", corOptions = "use = 'p'")
moduleOutput <- data.frame(names(datExpr))
moduleOutput[,2] <- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3] <- intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
write.table(moduleOutput, "VS_Net_ModuleOutput_grey.txt", sep="\t", quote=F)
moduleOutput1 <- moduleOutput[which(moduleOutput$ModuleColor != "grey"),]
write.table(moduleOutput1, "VS_Net_ModuleOutput_nogrey.txt", sep="\t", quote=F)

TOM <- TOMsimilarityFromExpr(datExpr, power= PWR, corType = "pearson", networkType = "signed", TOMType = "signed", TOMDenom = "mean", nThreads = 24, verbose = 5, indent = 0)
colnames(TOM) <- rownames(TOM) <- colnames(datExpr)
save(TOM, file="VS_Net_TOM_SIGNED.RData")
Connectivity <- apply(TOM, 1, sum)
save(Connectivity, file = "VS_Net_Connectivity.RData")

dir.create("Cyto")
setwd("Cyto")
for(module in unique(moduleColorsAutomatic))
    {
    inModule <- is.finite(match(moduleColorsAutomatic, module))
    modTOM <- TOM[inModule, inModule]
    cyt <- exportNetworkToCytoscape(modTOM, edgeFile=paste("CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile = paste("CytoNode", paste(module, collapse="-"), ".txt", sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColorsAutomatic[inModule], nodeNames = names(datExpr)[inModule])
    }
setwd("./../")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
df <- melt(MEs0)
dir.create("Plot")
setwd("Plot")
df$Rows <- factor(df$Rows, levels = colnames(tab))
df <- df[!df$variable == "MEgrey", ]

for (colr in unique(df$variable))
    {
    print(colr)
    PLOT<-ggplot(data=subset(df, variable == colr), aes(x=Rows, y=value)) + geom_bar(aes(fill=Class), stat="identity",position = "dodge") + scale_y_continuous(limits=c(-1,+1)) + theme_bw() + theme(strip.text.x = element_text(size=14, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF")) + scale_fill_manual(values = c("green", "blue", "red", "gold", "green4", "blue4", "red4", "gold4")) + theme(axis.title.x = element_blank(),axis.text.x  = element_text(face="bold", size=6,angle = 90, hjust = 1)) + theme(axis.title.y = element_blank(),axis.text.y  = element_text(face="bold", size=6)) + labs(title=colr)
    plotfile<-paste(colr, "barplot.pdf", sep="_")
    ggsave(plotfile, PLOT)
    }

setwd("./../")
save(df, MEs0, file="VS_Net_ModulesBarplotData.RData")


## module traits association heatmap
moduleColorsNEW <- moduleColorsAutomatic
MEs1 <- moduleEigengenes(datExpr, moduleColorsNEW)$eigengenes
MEsNEW <- MEs1
MEsNEW$MEgrey <-NULL


datTraits2 <- datTraits
datTraits2$HET <- as.numeric(gsub("2", 0, datTraits2$GENOTYPE))
datTraits2$WT <- as.numeric(gsub("1", 0, datTraits2$GENOTYPE))
datTraits2$WT <- as.numeric(gsub("2", 1, datTraits2$WT))

datTraits2$THREE <- as.numeric(gsub("2|3|4", 0, datTraits2$AGE))
datTraits2$SIX <- as.numeric(gsub("1|3|4", 0, datTraits2$AGE))
datTraits2$SIX <- as.numeric(gsub("2", 1, datTraits2$SIX))
datTraits2$TWELVE <- as.numeric(gsub("1|2|4", 0, datTraits2$AGE))
datTraits2$TWELVE <- as.numeric(gsub("3", 1, datTraits2$TWELVE))
datTraits2$NINETEEN <- as.numeric(gsub("1|2|3", 0, datTraits2$AGE))
datTraits2$NINETEEN <- as.numeric(gsub("4", 1, datTraits2$NINETEEN))

datTraits2$GENOTYPE <- NULL
datTraits2$AGE <- NULL

modTraitCor <- cor(MEsNEW, datTraits2, method="spearman") 
write.table(modTraitCor,"VS_Net_modTraitCor_New.txt", sep="\t", quote=F, row.names=T, col.names=T)

modTraitP <- corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\t(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) <- dim(modTraitCor) 

pdf("VS_Net_Heatmap_DatTraits_Pval.pdf", height = 18, width = 9)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = as.matrix(modTraitCor), xLabels = names(datTraits2), yLabels = names(MEsNEW), ySymbols = names(MEsNEW), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.65, zlim = c(-1,1), main = paste("Module Eigengene Correlation to Traits"))
dev.off()

modTraitCor <- as.data.frame(modTraitCor)
modTraitCorSorted <- modTraitCor[order(modTraitCor$HET, modTraitCor$WT, modTraitCor$THREE, modTraitCor$SIX, modTraitCor$TWELVE, modTraitCor$NINETEEN, decreasing = T),]
write.table(modTraitCorSorted,"VS_Net_modTraitCor_New_Sorted.txt", sep="\t", quote=F, row.names=T, col.names=T)
modTraitPSorted <- corPvalueStudent(as.matrix(modTraitCorSorted), nSamples)
textMatrixSorted = paste(signif(as.matrix(modTraitCorSorted), 2), "\t(",signif(modTraitPSorted, 1), ")", sep = "")
dim(textMatrixSorted) <- dim(modTraitCorSorted) 

pdf("VS_Net_Heatmap_DatTraits_Pval_sorted.pdf", height = 18, width = 12)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = as.matrix(modTraitCorSorted), xLabels = names(datTraits2), yLabels = row.names(modTraitCorSorted), ySymbols = row.names(modTraitCorSorted), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrixSorted, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1), main = paste("Module Eigengene Correlation to Traits"))
dev.off()

save.image(file =  "VS_WGCNA.RData")

```

&nbsp;

## Network Enrichment
```{R}

############################
## LOAD LIBRARIES
rm(list = ls())
library(SuperExactTest)
library(ggplot2)
library(ggnewscale)
library(rlist)
library(wesanderson)

############################
## LOAD MODULE OUTPUT
moduleOutput <- read.table("./../VS_Net_ModuleOutput_nogrey.txt", header = TRUE, sep = "\t")
moduleOutput$Module <- paste("ME", moduleOutput$ModuleColor, sep = "")

############################
## LOAD MODULE TRAIT CORR
moduleTrait <- read.table("./../VS_Net_modTraitCor_New.txt", header = TRUE, sep = "\t")
moduleTrait$Module <- row.names(moduleTrait)
moduleTrait <- moduleTrait[order(moduleTrait$WT, decreasing = TRUE),]
moduleTrait$NewMod <- paste("M", sprintf("%02d", seq(1:nrow(moduleTrait))), sep = "")
moduleTrait$NewModule <- paste(moduleTrait$NewMod, moduleTrait$Module, sep = " | ")

moduleOutputSortedTemp <- merge(moduleOutput, moduleTrait, by = "Module")
moduleOutputSorted <- moduleOutputSortedTemp[order(moduleOutputSortedTemp$NewModule),]


############################
## CREATE LIST OF MODULE GENES
moduleGenes <- vector("list", length(unique(sort(moduleOutputSorted$NewModule))))
names(moduleGenes) <- as.character(unique(sort(moduleOutputSorted$NewModule)))

myModData <- NULL

for(myMod in unique(sort(moduleOutputSorted$NewModule)))
  {
  myModOutput <- moduleOutputSorted[moduleOutputSorted$NewModule == myMod,]
  myModGenes <- as.character(unique(sort(myModOutput$Gene)))
  print(paste("Module: ", myMod, " | Genes: ", length(myModGenes), sep = "")) 
  moduleGenes[[myMod]] <- myModGenes

  myRow <- cbind(as.character(myMod), as.numeric(length(myModGenes)))
  myModData <- rbind(myModData, myRow)
  }

myModData <- data.frame(myModData)
colnames(myModData) <- c("Module", "Size")
row.names(myModData) <- myModData$Module

write.table(myModData, "VS_Net_ModuleSize_nogrey.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


############################
## CREATE SORTED LIST OF MODULE GENES
moduleGenesSorted <- vector("list", nrow(myModData))
names(moduleGenesSorted) <- myModData$Module

for(newMod in unique(sort(moduleOutputSorted$NewModule)))
  {
  # newModOutput <- moduleOutput[moduleOutput$ModuleColor == newMod,]
  newModOutput <- moduleOutputSorted[grepl(newMod, moduleOutputSorted$NewModule),]
  newModGenes <- as.character(unique(sort(newModOutput$Gene)))
  print(paste("Module: ", newMod, " | Genes: ", length(newModGenes), sep = "")) 
  moduleGenesSorted[[newMod]] <- newModGenes
  }



############################
## LOAD DEG LISTS
vs_three_up <- scan("VS_03_Months_DEG_SIG_GENES_UP.txt", what = "", sep = "\n")    #  49
vs_three_dn <- scan("VS_03_Months_DEG_SIG_GENES_DN.txt", what = "", sep = "\n")    #  58
vs_six_up <- scan("VS_06_Months_DEG_SIG_GENES_UP.txt", what = "", sep = "\n")      #  74
vs_six_dn <- scan("VS_06_Months_DEG_SIG_GENES_DN.txt", what = "", sep = "\n")      # 201
vs_twelve_up <- scan("VS_12_Months_DEG_SIG_GENES_UP.txt", what = "", sep = "\n")   #  34
vs_twelve_dn <- scan("VS_12_Months_DEG_SIG_GENES_DN.txt", what = "", sep = "\n")   # 107
vs_nineteen_up <- scan("VS_19_Months_DEG_SIG_GENES_UP.txt", what = "", sep = "\n") #  48
vs_nineteen_dn <- scan("VS_19_Months_DEG_SIG_GENES_DN.txt", what = "", sep = "\n") #  56


############################
## LOAD CHANNEL GENE LIST
channel_genes <- c("Kcnn2", "Kcnc3", "Kcnd3", "Kcnma1", "Itpr1", "Trpc3", "Orai2", "Stim1", "Slc8a1", "Slc8a2", "Slc8a3", "Cacna1g", "Cacna1a", "Atp2b2")



############################
## COMBINED GENES LIST
deGenes <- list("01_03_Months_UP" = vs_three_up, 
                "02_03_Months_DN" = vs_three_dn, 
                "03_06_Months_UP" = vs_six_up, 
                "04_06_Months_DN" = vs_six_dn, 
                "05_12_Months_UP" = vs_twelve_up, 
                "06_12_Months_DN" = vs_twelve_dn, 
                "07_19_Months_UP" = vs_nineteen_up, 
                "08_19_Months_DN" = vs_nineteen_dn, 
                "09_Chennel_Genes" = channel_genes)


############################
## PERFORM supertest
setout <- list()
i <- 0

for(mod in names(moduleGenesSorted))
 {
 i <- i+1
 newGenes <- list.append(deGenes, mod = moduleGenesSorted[[mod]])
 names(newGenes)[[length(newGenes)]] <- mod
 setres <- supertest(newGenes, n = 13000, degree = 2)
 setresData <- as.data.frame(summary(setres)$Table)
 setresDataSel <- setresData[grep(mod, setresData$Intersections),c(1, 3, 5, 6)]
 setout[[i]] <- setresDataSel
 print(paste(mod, i, sep = " "))
 }

names(setout) <- names(moduleGenesSorted)

setoutData <- Reduce("rbind", setout)

save(setout, setoutData, file = "VS_Net_SET.RData")


############################
## DATA FOR PLOTS
setoutDataTemp <- as.data.frame(matrix(unlist(strsplit(setoutData$Intersections, " & ")), byrow = T, ncol = 2))
colnames(setoutDataTemp) <- c("DEG_PAIR", "Module")
setoutDataTemp$Intersections <- setoutData$Intersections

setoutData2 <- merge(setoutData, setoutDataTemp, by = "Intersections")
setoutData2$degpair <- as.numeric(substr(setoutData2$DEG_PAIR, 1, 2)) + 2
setoutData2$NegLog10 <- -log10(setoutData2$P.value + 10^-10)

setoutData3 <- setoutData2

write.table(setoutData3, "VS_Net_SET_1.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

save(moduleOutput, moduleOutputSorted, moduleTrait, myModData, moduleGenesSorted, deGenes, setoutData3, file = "VS_Net_SET_DATA.RData")






############################
## LOAD LIBRARIES
rm(list = ls())
library(SuperExactTest)
library(ggplot2)
library(ggnewscale)
library(rlist)
library(wesanderson)
library(reshape2)

load("VS_Net_SET_DATA.RData")
# [1] "deGenes"            "moduleGenesSorted"  "moduleOutput"      
# [4] "moduleOutputSorted" "moduleTrait"        "myModData"         
# [7] "setoutData3"

moduleTraitCor <- melt(moduleTrait)
colnames(moduleTraitCor) <- c("Old_Module", "NewMod", "Module", "State", "Corr")
moduleTraitCor$New <- rep(0, nrow(moduleTraitCor))
moduleTraitCor[moduleTraitCor$State == "WT", "New"] <- 12
moduleTraitCor[moduleTraitCor$State == "HET", "New"] <- 13
moduleTraitCor[moduleTraitCor$State == "THREE", "New"] <- 14
moduleTraitCor[moduleTraitCor$State == "SIX", "New"] <- 15
moduleTraitCor[moduleTraitCor$State == "TWELVE", "New"] <- 16
moduleTraitCor[moduleTraitCor$State == "NINETEEN", "New"] <- 17

############################
## COLOR UP & DN SEPERATELY
setoutData32 <- setoutData3
setoutData32$PAIR <- rep("MISC", nrow(setoutData32))
setoutData32[grepl("^01", setoutData32$DEG_PAIR),"PAIR"] <- "01_03_Months_UP"
setoutData32[grepl("^02", setoutData32$DEG_PAIR),"PAIR"] <- "02_03_Months_DN"
setoutData32[grepl("^03", setoutData32$DEG_PAIR),"PAIR"] <- "03_06_Months_UP"
setoutData32[grepl("^04", setoutData32$DEG_PAIR),"PAIR"] <- "04_06_Months_DN"
setoutData32[grepl("^05", setoutData32$DEG_PAIR),"PAIR"] <- "05_12_Months_UP"
setoutData32[grepl("^06", setoutData32$DEG_PAIR),"PAIR"] <- "06_12_Months_DN"
setoutData32[grepl("^07", setoutData32$DEG_PAIR),"PAIR"] <- "07_19_Months_UP"
setoutData32[grepl("^08", setoutData32$DEG_PAIR),"PAIR"] <- "08_19_Months_DN"
setoutData32[grepl("^09", setoutData32$DEG_PAIR),"PAIR"] <- "09_Chennel_Genes"


setoutData32$NGL10 <- setoutData32$NegLog10
setoutData32[setoutData32$NegLog10 < -log10(0.05 + 10^-10), "NGL10"] <- 0
setoutData32[grepl("DN", setoutData32$DEG_PAIR), "NGL10"] <- -1 * setoutData32[grepl("DN", setoutData32$DEG_PAIR),"NGL10"]

# setoutData32$Module <- factor(setoutData32$Module, levels = myModDataSorted$Module)



############################
## CUSTOM FUNCTION TO ADD 2 COLOR SCALES

ggplot_add.new_aes <- function(object, plot, object_name) 
  {
  plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
  plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
  plot$labels <- bump_aes(plot$labels, new_aes = object)
  plot
  }

new_scale <- function(new_aes) 
  {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
  }






############################
## SELECT ALL MODULES
pAll <- ggplot(data = setoutData32) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "01_03_Months_UP",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        scale_fill_gradient2(low = "green2", mid = "grey95", high = "red2") +
        ylim(c(0, max(setoutData32[setoutData32$DEG_PAIR == "01_03_Months_UP","degpair"]) + 17)) + 
        # coord_polar(theta="x") +
        theme(panel.background = element_blank(), # bg of the panel
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_text(size=10, angle=90, vjust = 0.5, hjust=1),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=0)) +
        NULL

pAll <- pAll + 
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "02_03_Months_DN",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "03_06_Months_UP",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "04_06_Months_DN",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "05_12_Months_UP",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "06_12_Months_DN",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "07_19_Months_UP",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "08_19_Months_DN",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        ## Channel Genes
        new_scale("fill") + 
        geom_tile(data = setoutData32[setoutData32$DEG_PAIR == "09_Chennel_Genes",], aes(x = Module, y = degpair, fill = NGL10), color = "white", size = 0.5) +
        scale_fill_gradient(low = "grey95", high = "blue") +
        ## MODULE TRAIT CORRELATION
        new_scale("fill") + 
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "WT",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        scale_fill_gradient2(low = "blue", mid = "white", high = "red2") +
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "HET",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "THREE",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "SIX",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "TWELVE",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        geom_tile(data = moduleTraitCor[moduleTraitCor$State == "NINETEEN",], aes(x = Module, y = New, fill = Corr), color = "white", size = 0.5) + 
        ## ADDING LINES TO SEPERATE
        geom_hline(yintercept = 2.5) + 
        geom_hline(yintercept = 4.5) + 
        geom_hline(yintercept = 6.5) + 
        geom_hline(yintercept = 8.5) + 
        geom_hline(yintercept = 10.5) +
        geom_hline(yintercept = 11.5) +
        geom_hline(yintercept = 13.5) +
        geom_hline(yintercept = 17.5) +
        NULL

ggsave(filename = "VS_Net_SET_1.pdf", plot = pAll, width = 16, height = 12, dpi = 300)


pAll <- pAll + theme(legend.position = "none")
ggsave(filename = "VS_Net_SET_2.pdf", plot = pAll, width = 16, height = 12, dpi = 300)

```

&nbsp;


## Plot Network
```{R}
rm(list = ls())
library(Seurat)
library(reshape2)
library(future)
library(furrr)
library(purrr)
library(ComplexHeatmap)
library(scales)
library(cowplot)
library(circlize)
library(dplyr)
library(pheatmap)
library(igraph)
library(rgexf)
library(qgraph)

layoutlist <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

plotNet <- function(ifn, top, prefx, colr)
    {
    # Read Data
    print(paste("Processing", prefx, "for Top", top, "connections", sep = " "))
    module <- read.table(ifn, sep = "\t", header = TRUE)
    moduleSorted <- module[order(module$weight, decreasing = TRUE),]
    moduleSortedTop <- moduleSorted[1:top,]
    moduleSortedTop$direction <- NULL
    moduleSortedTop$fromAltName <- NULL
    moduleSortedTop$toAltName <- NULL
    colnames(moduleSortedTop) <- c("Source", "Target", "Weight")

    # Nodes
    nodes <- as.data.frame(unique(sort(c(as.character(moduleSortedTop$Source), as.character(moduleSortedTop$Target)))))
    colnames(nodes) <- "id"

    # Edges
    edges <- moduleSortedTop
    colnames(edges) <- c("from", "to", "weight")
    write.table(edges, paste(prefx, "_Top_", top, "_Edges.tsv", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")

    # Plot Network
    g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    e <- get.edgelist(g, names = FALSE)
    d <- degree(g)
    l <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(g))
    plot(g, layout = l, vertex.size = d, vertex.label.cex = 0.5, vertex.label.family = "Helvetica")

    colrs <- c("gray50", "tomato", "gold")

    V(g)$color <- colr
    V(g)$frame.color <- "white"
    V(g)$size <- d*1.2
    V(g)$label.family <- "Helvetica"
    V(g)$label.color <- "black"
    V(g)$label.cex <- 0.75
    V(g)$label.font <- 3
    E(g)$edge.color <- "gray80"
    E(g)$curved <- 0
    # graph_attr(g, "layout") <- layout_on_sphere
    l <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), area = 8*(vcount(g)^2), repulse.rad = (vcount(g)^2.5))
    pdf(paste(prefx, "_Top_", top, "_1.pdf", sep = ""), width = 9, height = 9)
    plot(g, layout = l)
    dev.off()


    E(g)$curved <- 0.5
    pdf(paste(prefx, "_Top_", top, "_2.pdf", sep = ""), width = 9, height = 9)
    plot(g, layout = l)
    dev.off()


    V(g)$label.cex <- d/18
    V(g)$label.color <- "white"
    E(g)$curved <- 0
    pdf(paste(prefx, "_Top_", top, "_3.pdf", sep = ""), width = 9, height = 9)
    plot(g, layout = l)
    dev.off()

    }


## Module: darkmagenta
top <- 250
ifn <- "/work/psychiatry/rgree2/BULKSEQ/L_WGCNA/G_WGCNA_02/Cyto/CytoEdgedarkmagenta.txt"
prefx <- "Module_DarkMagenta"
colr <- "darkmagenta"
plotNet(ifn, top, prefx, colr)


## Module: brown
top <- 250
ifn <- "/work/psychiatry/rgree2/BULKSEQ/L_WGCNA/G_WGCNA_02/Cyto/CytoEdgebrown.txt"
prefx <- "Module_Brown"
colr <- "brown"
plotNet(ifn, top, prefx, colr)



## Module: darkred
top <- 250
ifn <- "/work/psychiatry/rgree2/BULKSEQ/L_WGCNA/G_WGCNA_02/Cyto/CytoEdgedarkred.txt"
prefx <- "Module_DarkRed"
colr <- "darkred"
plotNet(ifn, top, prefx, colr)



## Module: skyblue
top <- 250
ifn <- "/work/psychiatry/rgree2/BULKSEQ/L_WGCNA/G_WGCNA_02/Cyto/CytoEdgeskyblue.txt"
prefx <- "Module_SkyBlue"
colr <- "skyblue"
plotNet(ifn, top, prefx, colr)


```