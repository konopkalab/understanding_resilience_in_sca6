# DIFFERENTIAL GENE EXPRESSION

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
