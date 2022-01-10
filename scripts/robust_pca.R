
# ROBUST PCA ANALYSIS -----------------------------------------------------

# Here we aim to identify outlier samples that are not obviously seen in common PCA.
# We use the methodology implemented in the 'rrcov' package.

library(rrcov)
library(DESeq2)

# Load count table and annotation data
load("results/txi/txi_gene.rda")
load("results/important_variables/ann.rda")

# DESeq2 object
dds <- DESeqDataSetFromTximport(txi,
                                colData = ann,
                                design = ~ group)

v <- vst(dds, blind = F)
v <- t(assay(v))
pc1 <- PcaGrid(v, 5)

# Plot rPCA
pdf("results/plots_paper/robust_pca.pdf")
plot(pc1)
dev.off()

# The SRR5961961 sample is identified as divergent from others, and for this reason is excluded from 
# downstream analyses. 

# plot(pc1$scores[,1], pc1$od)
# text(pc1$scores[pc1$od > 80, 1], pc1$od[pc1$od > 80]+2, labels = rownames(pc1$scores)[which(pc1$od > 80)])
# pc2 <- PcaHubert(v, 5)
# plot(pc2$scores[,1], pc2$od)
# text(pc2$scores[pc2$od > 90, 1], pc2$od[pc2$od > 90]-2, labels = rownames(pc2$scores)[which(pc2$od > 90)])
