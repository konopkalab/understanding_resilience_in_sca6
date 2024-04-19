# DIFFERENTIAL GENE EXPRESSION ANALYSIS

## Combine count tables, generate CPM and RPKM tables
```{R}
## load libraries
rm(list = ls())
library(GenomicFeatures)
library(rtracklayer)
library(edgeR)

## load count dataset
files <- list.files(path = "path_to_counts_directory", pattern = ".txt")
samples <- gsub("_OUTPUT.LastLinesRem.txt", "", files )
myfiles <- lapply(paste("path_to_counts_directory", files, sep = ""), read.table, sep = "\t", fill = TRUE)
count <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "V1") }, myfiles)
rownames(count) <- count$V1
count$V1 <- NULL
names(count) <- samples

## load gene lengths
df <- read.table("gencode.vM17.protein_coding_length.txt", header = T, row.names = 1)

## merge gene lengths with combined count tables
new <- merge(df, count, by = "row.names")
geneLengths <- as.data.frame(new[,2])
row.names(geneLengths) <- new[,1]
names(geneLengths) <- "Length"

## filter and save raw counts table
countTemp <- count[row.names(count) %in% row.names(geneLengths),]
counts <- countTemp
write.table(counts, "VS_COUNT_TABLE.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## calculate and save RPKM table
rpkm <- rpkm(counts, log = FALSE, gene.length = geneLengths$Length)
write.table(rpkm, "VS_RPKM_TABLE.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## calculate and save CPM table
cpm <- cpm(counts, log = FALSE)
write.table(cpm, "VS_CPM_TABLE.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## load meta data table
metaTempA <- read.table("VS_META.txt", header = TRUE, sep = "\t")
row.names(metaTempA) <- metaTempA$Sample

metaTempB <- as.data.frame(matrix(unlist(strsplit(colnames(counts), "_")), ncol = 4, byrow = TRUE))
row.names(metaTempB) <- colnames(counts)
colnames(metaTempB) <- c("Project", "Animal", "Replicate", "SeqInfo")

## update meta data
meta <- merge(metaTempA, metaTempB, by = "row.names")
row.names(meta) <- meta$Row.names
meta$Row.names <- NULL

## write updated meta data table
write.table(meta, "VS_META_TABLE.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## save counts, cpm, rpkm and meta data as RData
save(counts, rpkm, cpm, meta, file = "VS_COUNTS_RPKM_CPM_META.RData")

```


&nbsp;

## Principal Component Analysis and Hierarchical Clustering
```{R}
## load libraries
rm(list = ls())
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(edgeR)
library(ggrepel)
library(ggfortify)
library(dendextend)

## read counts, cpm, rpkm and meta data
listAll <- load("VS_COUNTS_RPKM_CPM_META.RData")

## read updated meta data
newmeta <- read.table("VS_META.txt", sep = "\t", header = TRUE)
row.names(newmeta) <- newmeta$SAMPLE_1

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


## subset for specific age and genotype
myage <- c("03m")
META <- as.data.frame(newmeta[newmeta$AGE %in% myage,])
COUNTS <- as.data.frame(t(counts2[row.names(counts2) %in% row.names(META),]))
RPKM <- as.data.frame(t(rpkm2[row.names(rpkm2) %in% row.names(META),]))
CPM <- as.data.frame(t(cpm2[row.names(cpm2) %in% row.names(META),]))


## PCA for RPKM
pca1 <- prcomp(t(log2(RPKM + 1)))
pca1DataTemp <- as.data.frame(pca1$x)
pca1Data <- merge(pca1DataTemp, META, by = "row.names", all = TRUE)

pca1Data$SAMPLE <- as.factor(pca1Data$SAMPLE)
pca1Data$GENOTYPE <- as.factor(pca1Data$GENOTYPE)
pca1Data$AGE <- as.factor(pca1Data$AGE)

pcaPlot1 <- ggplot(pca1Data, aes(PC1, PC2)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            # geom_text_repel(aes(PC1, PC2, label = Sample)) + 
            theme_bw() + 
            labs(title = "PCA | log2 RPKM") + 
            NULL
ggsave(filename = "VS_RPKM_PCA_PC1_PC2_Genotype.pdf", plot = pcaPlot1, width = 7, height = 6, dpi = 300)  

pcaPlot1 <- ggplot(pca1Data, aes(PC1, PC3)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            # geom_text_repel(aes(PC1, PC2, label = Sample)) + 
            theme_bw() + 
            labs(title = "PCA | log2 RPKM") + 
            NULL
ggsave(filename = "VS_RPKM_PCA_PC1_PC3_Genotype.pdf", plot = pcaPlot1, width = 7, height = 6, dpi = 300)  

pcaPlot1 <- ggplot(pca1Data, aes(PC2, PC3)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            # geom_text_repel(aes(PC1, PC2, label = Sample)) + 
            theme_bw() + 
            labs(title = "PCA | log2 RPKM") + 
            NULL
ggsave(filename = "VS_RPKM_PCA_PC2_PC3_Genotype.pdf", plot = pcaPlot1, width = 7, height = 6, dpi = 300)  


## PCA for CPM
pca2 <- prcomp(t(log2(CPM + 1)))
pca2DataTemp <- as.data.frame(pca2$x)
pca2Data <- merge(pca2DataTemp, META, by = "row.names", all = TRUE)

pca2Data$SAMPLE <- as.factor(pca2Data$SAMPLE)
pca2Data$GENOTYPE <- as.factor(pca2Data$GENOTYPE)
pca2Data$AGE <- as.factor(pca2Data$AGE)

pcaPlot2 <- ggplot(pca2Data, aes(PC1, PC2)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            theme_bw() + 
            labs(title = "PCA | log2 CPM") + 
            NULL
ggsave(filename = "VS_CPM_PCA_PC1_PC2_Genotype.pdf", plot = pcaPlot2, width = 7, height = 6, dpi = 300)  

pcaPlot2 <- ggplot(pca2Data, aes(PC1, PC3)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            theme_bw() + 
            labs(title = "PCA | log2 CPM") + 
            NULL
ggsave(filename = "VS_CPM_PCA_PC1_PC3_Genotype.pdf", plot = pcaPlot2, width = 7, height = 6, dpi = 300)  

pcaPlot2 <- ggplot(pca2Data, aes(PC2, PC3)) + 
            geom_point(aes(shape = GENOTYPE, colour = GENOTYPE), size = 3) + 
            stat_ellipse(geom="polygon", aes(fill = GENOTYPE), alpha = 0.2, show.legend = FALSE, level = 0.95) +
            guides(
              colour = guide_legend(order = 1), 
              shape = guide_legend(order = 2)) +
            theme_bw() + 
            labs(title = "PCA | log2 CPM") + 
            NULL
ggsave(filename = "VS_CPM_PCA_PC2_PC3_Genotype.pdf", plot = pcaPlot2, width = 7, height = 6, dpi = 300)  


## Save RData
save(pca1, pca1Data, pcaPlot1, pca2, pca2Data, pcaPlot2, file = "VS_PCA.RData")


## HCLUST for RPKM
RPKMT <- as.data.frame(t(log2(RPKM+1)))

dendRPKM <- RPKMT %>% scale %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k = 2) %>% set("branches_lwd", 1.5)
ggd1 <- as.ggdend(dendRPKM)
plotDendro1 <- ggplot(ggd1, horiz = TRUE) + theme(plot.margin=unit(c(1,5,1,1),"cm"))
ggsave(filename = "VS_RPKM_HCLUST.pdf", plot = plotDendro1, width = 16, height = 8, units = "in", dpi = 300)


## HCLUST for CPM
CPMT <- as.data.frame(t(log2(CPM+1)))

dendCPMT <- CPMT %>% scale %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k = 2) %>% set("branches_lwd", 1.5)
ggd2 <- as.ggdend(dendCPMT)
plotDendro2 <- ggplot(ggd2, horiz = TRUE) + theme(plot.margin=unit(c(1,5,1,1),"cm"))
ggsave(filename = "VS_CPM_HCLUST.pdf", plot = plotDendro2, width = 16, height = 8, units = "in", dpi = 300)


## Save RData
save(RPKMT, dendRPKM, plotDendro1, CPMT, dendCPMT, plotDendro2, file = "VS_HCLUST.RData")


setwd("./../")

## Repeat for other ages (6m, 12m, 19m)
## Repeat for both the genotypes (wt and het)

```

&nbsp;

## Differential Gene Expression using DESeq2
```{R}
## load libraries
rm(list = ls())
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

## mito genes
genesMito <- scan("mtGenes.txt", what = "", sep = "\n")

## read counts, cpm, rpkm and meta data
listAll <- load("VS_COUNTS_RPKM_CPM_META.RData")

## load updated meta data
newmeta <- read.table("/work/Neuroinformatics_Core/akulk1/FOR_NEUROINFO_CORE/VIKRAM/VS_META.txt", sep = "\t", header = TRUE)
row.names(newmeta) <- newmeta$SAMPLE_1

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

DESIGN <- newmeta
COUNTS <- as.data.frame(t(counts2[row.names(counts2) %in% row.names(DESIGN),!colnames(counts2) %in% genesMito]))
RPKM <- as.data.frame(t(rpkm2[row.names(rpkm2) %in% row.names(DESIGN),!colnames(rpkm2) %in% genesMito]))
CPM <- as.data.frame(t(cpm2[row.names(cpm2) %in% row.names(DESIGN),!colnames(cpm2) %in% genesMito]))

## Filtering parameters
cutoff_l2fc <- 0.3 ## log2(1.2) = 0.26 | log2(1.25) = 0.3
cutoff_p <- 0.05

## DGE parameters
cond1 <- "WT"
cond2 <- "HET"
type <- "RPKM"
age <- "03m"

## filter outliers samples and genes with low expression and run DESeq2
prefix <- paste(cond1, "_X_", cond2, "_", age, "_", type, sep = "")
print(prefix)

outliers <- c("WT_03m_SCA6_35_2", "HET_03m_SCA6_34_6")

dfSelectedTemp <- RPKM[,!colnames(RPKM) %in% outliers]
designSelectedTemp <- DESIGN[!row.names(DESIGN) %in% outliers,]

designSelectedTemp2 <- designSelectedTemp[designSelectedTemp$GENOTYPE %in% c(cond1, cond2),]
designSelected <- designSelectedTemp[designSelectedTemp$AGE %in% c(age),]

dfSelected <- dfSelectedTemp[,colnames(dfSelectedTemp) %in% row.names(designSelected)]
dfSelected <- dfSelected[,order(colnames(dfSelected))]

filter <- apply(dfSelected, 1, function(x) { all(x[1:4] >= 0.5) | all(x[5:8] >= 0.5)})
dfFiltered <- dfSelected[filter,]
myCountsTemp <- COUNTS[filter, colnames(COUNTS) %in% row.names(designSelected)]
myCounts <- myCountsTemp[,order(colnames(myCountsTemp))]
myDesign <- designSelected[order(row.names(designSelected)),]
print(dim(myCounts))

ddsNew <- DESeqDataSetFromMatrix(myCounts, colData = myDesign, design = ~ SEX + RIN + GENOTYPE)
ddsNew <- DESeq(ddsNew, test = "Wald", fitType = "mean")
resNew <- results(ddsNew, contrast=c("GENOTYPE", cond2, cond1))
resNew2 <- as.data.frame(resNew)
resNew2 <- na.omit(resNew2)

resNew3 <- resNew2[resNew2$padj <= cutoff_p,]
resNew4 <- resNew3[resNew3$log2FoldChange >= cutoff_l2fc | resNew3$log2FoldChange <= -(cutoff_l2fc),]
resNew4up <- resNew4[resNew4$log2FoldChange >= cutoff_l2fc,]
resNew4dn <- resNew4[resNew4$log2FoldChange < -cutoff_l2fc,]

print(dim(resNew4))
print(dim(resNew4up))
print(dim(resNew4dn))

## heatmap for significantly differentially expressed genes
degRPKM <- dfFiltered[row.names(resNew4),]
pheatmap(degRPKM, scale="row", cluster_rows=T, cluster_cols=T, border_color=NA, show_rownames=F, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), filename = paste("VS_DEG_", prefix, "_Heatmap.pdf", sep = ""), width = 6, height = 6)

## MA plot
MAplot <- ggplot(data = resNew2, aes(x = log2(baseMean+1), y = log2FoldChange)) +
          geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
          geom_point(data=resNew4up, colour="green3", size=1.5, shape=16, alpha=0.8) +
          geom_point(data=resNew4dn, colour="royalblue", size=1.5, shape=16, alpha=0.8) +
          annotate(geom = "text", x = 4, y = 2, label = nrow(resNew4up)) +
          annotate(geom = "text", x = 4, y = -2, label = nrow(resNew4dn)) +
          theme_bw() +
          ylim(-2.5,2.5) +
          labs(x="Log2 Mean Expression", y="Log2 Fold Change") +
          NULL
ggsave(filename = paste("VS_DEG_", prefix, "_MAPlot.pdf", sep = ""), plot = MAplot, width = 6, height = 6, units = "in", dpi = 300, useDingbats = FALSE)

### volcano plot
resNew2$neglogpval <- -log10(resNew2$padj + 10^-10)
resNew4up$neglogpval <- -log10(resNew4up$padj + 10^-10)
resNew4dn$neglogpval <- -log10(resNew4dn$padj + 10^-10)

VOLplot <- ggplot(data = resNew2, aes(x = log2FoldChange, y = neglogpval)) +
           geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
           geom_point(data=resNew4up, colour="green3", size=1.5, shape=16, alpha=0.8) +
           geom_point(data=resNew4dn, colour="royalblue", size=1.5, shape=16, alpha=0.8) +
           annotate(geom = "text", x = 2, y = 8, label = nrow(resNew4up)) +
           annotate(geom = "text", x = -2, y = 8, label = nrow(resNew4dn)) +
           theme_bw() +
           geom_vline(xintercept = log2(1.20), colour="red", linetype=2) + geom_vline(xintercept = -log2(1.20), colour="red", linetype=2) +
           geom_hline(yintercept = -log10(0.05), colour="red", linetype=2) +
           xlim(-2.5,2.5) +
           labs(x="log2 Fold Change", y="-log10 P-value")
ggsave(filename = paste("VS_DEG_", prefix, "_VolcanoPlot.pdf", sep = ""), plot = VOLplot, width = 6, height = 6, units = "in", dpi = 300, useDingbats = FALSE)

## write tables, lists and save output
write.table(resNew2, paste("VS_DEG_", prefix, ".txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(resNew4, paste("VS_DEG_", prefix, "_SIG.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(row.names(resNew4up), paste("VS_DEG_", prefix, "_SIG_GENES_UP.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(row.names(resNew4dn), paste("VS_DEG_", prefix, "_SIG_GENES_DN.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(degRPKM, paste("VS_DEG_", prefix, "_SIG_RPKM.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")
save(COUNTS, RPKM, DESIGN, filter, myCounts, myDesign, ddsNew, resNew, file = paste("VS_DEG_", prefix, ".RData", sep = ""))


## Repeat for other ages (6m, 12m, 19m)

```
