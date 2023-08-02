if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

counts <- as.matrix(read.csv("gencode_srr.csv", header = TRUE, row.names=1, sep ='\t'))
colData <- read.csv("sras_types.csv", sep = '\t', header = FALSE,row.names=1)
colData_breast <- colData
colData_blood <- colData
colData_brain <- colData
colData_breast[colData_breast == "blood"] <- "not breast"
colData_breast[colData_breast == "brain"] <- "not breast"

colData_brain[colData_brain == "blood"] <- "not brain"
colData_brain[colData_brain == "breast"] <- "not brain"

colData_blood[colData_blood == "brain"] <- "not blood"
colData_blood[colData_blood == "breast"] <- "not blood"

dds_breast <- DESeqDataSetFromMatrix(counts, colData = colData_breast, design = ~V2)
dds_breast <- DESeq(dds_breast)
res_breast <- results(dds_breast, alpha=0.05, lfcThreshold = log2(4), altHypothesis = "less")
res_breast_significant <- subset(res_breast, res_breast$padj < 0.05)
write.table(row.names(res_breast_significant), "breast_overexp.txt", sep = "\t", row.names = F, quote = F, col.names = F)

dds_brain <- DESeqDataSetFromMatrix(counts, colData = colData_brain, design = ~V2)
dds_brain <- DESeq(dds_brain)
res_brain <- results(dds_brain, alpha=0.05, lfcThreshold = log2(4), altHypothesis = "less")
res_brain_significant <- subset(res_brain, res_brain$padj < 0.05)
write.table(row.names(res_brain_significant), "brain_overexp.txt", sep = "\t", row.names = F, quote = F, col.names = F)

dds_blood <- DESeqDataSetFromMatrix(counts, colData = colData_blood, design = ~V2)
dds_blood <- DESeq(dds_blood)
res_blood <- results(dds_blood, alpha=0.05, lfcThreshold = log2(4), altHypothesis = "less")
res_blood_significant <- subset(res_blood, res_blood$padj < 0.05)
write.table(row.names(res_blood_significant), "blood_overexp.txt", sep = "\t", row.names = F, quote = F, col.names = F)

