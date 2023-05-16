if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("pheatmap")
library("pheatmap")


counts <- as.matrix(read.csv("deseq_breast.csv", header = TRUE, row.names=1, sep ='\t'))
colData <- read.csv("C:true_breast.csv", sep = ';', row.names=1)
colData$condition <- factor(colData$condition)

dds <- DESeqDataSetFromMatrix(counts, colData = colData, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
resOrderd <- res[order(res$padj),]
number_de_genes = sum(resOrderd$padj < 0.05, na.rm=TRUE)
res_de <- resOrderd[1:number_de_genes,]
res_de <- subset(res_de, abs(res_de$log2FoldChange) >= 2)

rescounts23 <- as.matrix(read.csv("deseq_breast_23.csv", header = TRUE, row.names=1, sep ='\t'))
dds23 <- DESeqDataSetFromMatrix(counts23, colData = colData, design = ~condition)
dds23 <- DESeq(dds23)
res23 <- results(dds23, alpha=0.05)
resOrderd23 <- res23[order(res23$padj),]
number_de_genes23 = sum(resOrderd23$padj < 0.05, na.rm=TRUE)
res_de23 <- resOrderd[1:number_de_genes23,]
res_de23 <- subset(res_de23, abs(res_de23$log2FoldChange) >= 2)

counts39 <- as.matrix(read.csv("deseq_breast_39.csv", header = TRUE, row.names=1, sep ='\t'))
dds39 <- DESeqDataSetFromMatrix(counts39, colData = colData, design = ~condition)
dds39 <- DESeq(dds39)
res39 <- results(dds39, alpha=0.05)
resOrderd39 <- res39[order(res39$padj),]
number_de_genes39 = sum(resOrderd39$padj < 0.05, na.rm=TRUE)
res_de39 <- resOrderd[1:number_de_genes39,]
res_de39 <- subset(res_de39, abs(res_de39$log2FoldChange) >= 2)

# Differentially expressed genes according to https://doi.org/10.1016/j.dib.2018.03.079
# genes = genes_up + genes_down
genes_up <- c("ADAMDEC1","GBP5", "IDO1", "SLAMF7", "CD274", "STAT1", "UBD", "MS4A1", "BANK1", "PLA2G2D", "CD8A", "ITGAL", "CD8B", "TIGIT", "IKZF3", "CD3G", "CXCL13", "FYB1", "SELL", "SLAMF6", "SH2D1A", "CYTIP", "IL7R", "ITK", "IKZF1", "PTPRC", "ITGA4", "CD3D", "CD2", "KLRK1", "CR2", "FGF10", "MMP9", "CD79A", "TBC1D10C", "CD19", "CD5", "CCR7", "LCK", "ITGB7", "HSH2D", "LILRB4", "SIT1", "CD27", "SPOCK2", "CD6", "CD247", "CARD11", "CD79B", "PRF1")
genes_down <- c("HMGCS2", "G0S2", "LIPE", "CYP4F22", "DGAT2", "ANGPTL4", "CD36", "AKR1C1", "FABP4", "AKR1C2", "ADIPOQ", "PLIN1", "LEP", "GPD1", "PLIN4", "AQP7", "CIDEC")
genes <- c("ADAMDEC1","GBP5", "IDO1", "SLAMF7", "CD274", "STAT1", "UBD", "MS4A1", "BANK1", "PLA2G2D", "CD8A", "ITGAL", "CD8B", "TIGIT", "IKZF3", "CD3G", "CXCL13", "FYB1", "SELL", "SLAMF6", "SH2D1A", "CYTIP", "IL7R", "ITK", "IKZF1", "PTPRC", "ITGA4", "CD3D", "CD2", "KLRK1", "CR2", "FGF10", "MMP9", "CD79A", "TBC1D10C", "CD19", "CD5", "CCR7", "LCK", "ITGB7", "HSH2D", "LILRB4", "SIT1", "CD27", "SPOCK2", "CD6", "CD247", "CARD11", "CD79B", "PRF1","HMGCS2", "G0S2", "LIPE", "CYP4F22", "DGAT2", "ANGPTL4", "CD36", "AKR1C1", "FABP4", "AKR1C2", "ADIPOQ", "PLIN1", "LEP", "GPD1", "PLIN4", "AQP7", "CIDEC")

# Cancer signature according to https://doi.org/10.3389/fgene.2022.912125
gene_sig <- c("KIF4A", "COL11A1", "RBP4", "SAA1", "SFRP1")

counts_map <- counts
colnames(counts_map) <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "N0", "N1", "N2")
pheatmap(t(counts_map[genes,]), cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")
