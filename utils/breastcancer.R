if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("pheatmap")
library("pheatmap")

library(cluster)
library(factoextra)

counts <- as.matrix(read.csv("Cutoff_1_Compressed_20_24_gene_expression_all.csv", header = TRUE, row.names=1, sep ='\t'))
colData <- read.csv("experiments_biosample.csv", sep = '\t', header = FALSE,row.names=1)

dds_all <- DESeqDataSetFromMatrix(counts, colData = colData, design = ~V2)
dds_all <- DESeq(dds_all)
res_all <- results(dds_all, alpha=0.05)
res_all_TN <- results(dds_all, alpha = 0.05, contrast=c("V2", "Triple Negative Breast Cancer Primary Tumor", "ER+ Breast Cancer Primary Tumor"),)

# Normalized counts
dds <- estimateSizeFactors(dds_all)
normalized_counts <- counts(dds_all, normalized=TRUE)

# TNBC analysis

# signature from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5469719/
# Not found:
genes_tnp <- c("PRKX", "UGT8", "HMGA1", "LPIN1", "HAPLN3", "FAM171A1", "BCL11A", "FOXC1", "ANKRD11")
# Not found: 
genes_tnp_under <- c("AR", "ANXA9","ATP8B1", "BCAS1", "CERS6","CCDC86", "CCDC96", "CMBL", "CXXC5", "DNAJC12", "ERBB4", "FAM174B", "GATA3", "LMNTD2-AS1", 
                     "MACIR","MLPH", "MMEL1", "NAT1", "TTC6", "PRR15", "RAB26", "RARA", "RHOB", "SLC44A4","SPATA20", "STARD10", "TBC1D9", "TFF1", "TFF3","TTC8", "TEC")

genes_tnp_all <- c(genes_tnp, genes_tnp_under)

# Look at DE analysis result of gene signature
# 2 where padj > 0.05, log2FoldChange is correct for all others
res_all_TN[genes_tnp_all,]

# Get Z-Scores
vsd <- assay(vst(dds_all))
Z <- t(scale(t(vsd)))
Z_tnp <- Z[genes_tnp_all, c(row.names(subset(colData, colData$V2 =="Triple Negative Breast Cancer Primary Tumor")),row.names(subset(colData, colData$V2 =="ER+ Breast Cancer Primary Tumor")))]
colnames(Z_tnp) <- c(paste0(rep("TNBC", 42), 1:42), paste0(rep("ERBC", 42), 1:42))
svg(filename="heatmap_tnbc.svg", width = 15, height = 3)
pheatmap(Z_tnp, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


# Clustering
colnames(Z_tnp) <- c(1:42, LETTERS, letters[1:16])
pam_res <- pam(t(Z_tnp), 2) # 4 misclassified

png(file="cluster_tnbc.png",   width     = 5.35,
    height    = 5.35,
    units     = "in",
    res       = 1000,
    pointsize = 2)
p<- fviz_cluster(pam_res, palette = "jco", legend.title = "Cluster", show.legend=F, main = "")
p <- p + geom_point(data=p$data[c(31,39),], aes(x=x, y=y), color=c("#0073C2"), shape = 2)
p + geom_point(data=p$data["m",], aes(x=x, y=y), color=c("#EFC000"), shape = 1)
dev.off()

#Signature from: https://www.nature.com/articles/s41467-017-01027-z
genes_tsa <- c("EGR1","EGR2", "FOS", "FOSB", "EGR3", "RCAN1", "TPPP", "NR4A3", "DPT", "CSRP1","JUND", "ATF3", "ACE",  "CXCL12", "CNN1", "FGL2", "MYADM", "CCN1")

res_all_nat <- results(dds_all, alpha = 0.05, contrast=c("V2", "Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor", "Reduction Mammoplasty - No known cancer"),)
res_er_nat <- results(dds_all, alpha = 0.05, contrast=c("V2", "ER+ Breast Cancer Primary Tumor", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor"),)
res_tn_nat <- results(dds_all, alpha = 0.05, contrast=c("V2", "Triple Negative Breast Cancer Primary Tumor", "Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor"),)
res_er_nocancer <- results(dds_all, alpha = 0.05, contrast=c("V2", "ER+ Breast Cancer Primary Tumor", "Reduction Mammoplasty - No known cancer"),)
res_nater_nocancer <-results(dds_all, alpha = 0.05, contrast=c("V2", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor", "Reduction Mammoplasty - No known cancer"),)
res_nat_nat <- results(dds_all, alpha = 0.05, contrast=c("V2", "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor", "Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor"),)

normalized_counts_tsa_er <- Z[genes_tsa,c(row.names(subset(colData, colData$V2 =="ER+ Breast Cancer Primary Tumor")),row.names(subset(colData, colData$V2 =="Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor")))]
colnames(normalized_counts_tsa_er) <- c(paste0(rep("", 42), 1:42), LETTERS, letters[1:4])

type <- c(rep("blue",42), rep("yellow",30))
pam_tsa_er <- pam(t(normalized_counts_tsa_er), 2) 
png(file="cluster_tsa_er.png",   width     = 5.35,
    height    = 5.35,
    units     = "in",
    res       = 1000,
    pointsize = 2)
p<- fviz_cluster(pam_tsa_er, palette = "jco", legend.title = "Cluster", show.legend=F, main = "")
p <- p + geom_point(data=p$data[c(2,24),], aes(x=x, y=y), color=c("#0073C2"), shape = 2)
p + geom_point(data=p$data[c("F","J","K","P","V","c"),], aes(x=x, y=y), color=c("#EFC000"), shape = 1)
dev.off()

normalized_counts_tsa_tn <- Z[genes_tsa,c(row.names(subset(colData, colData$V2 =="Triple Negative Breast Cancer Primary Tumor")),row.names(subset(colData, colData$V2 =="Uninvolved Breast Tissue Adjacent to TNBC Primary Tumor")))]
colnames(normalized_counts_tsa_tn) <- c(paste0(rep("", 42), 1:42), LETTERS[1:21])

pam_tsa_tn <- pam(t(normalized_counts_tsa_tn), 2) 
png(file="cluster_tsa_tn.png",   width     = 5.35,
    height    = 5.35,
    units     = "in",
    res       = 1000,
    pointsize = 2)
p<- fviz_cluster(pam_tsa_tn, palette = "jco", legend.title = "Cluster", show.legend=F, main = "")
p + geom_point(data=p$data[c(25,40,41),], aes(x=x, y=y), color=c("#0073C2"), shape = 2)
dev.off()
