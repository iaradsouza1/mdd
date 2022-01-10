
# Samples removal after robust pca analysis -------------------------------

library(dplyr)

# Load annotation 
load("results/important_variables/ann.rda")

outliers <- c("SRR5961961", "SRR5961809")

# Remove the 2 outliers samples from robust pca
ann %>% 
  tibble::rownames_to_column("run") %>% 
  filter(!(run %in% outliers)) -> ann
save(ann, file = "results/important_variables/ann.rda")

# Load txi information: gene
load("results/txi/txi_gene.rda")
txi$abundance <- txi$abundance[, !(colnames(txi$abundance) %in% outliers)]
txi$counts <- txi$counts[, !(colnames(txi$counts) %in% outliers)]
txi$length <- txi$length[, !(colnames(txi$length) %in% outliers)]
save(txi, file = "results/txi/txi_gene.rda")

# Load txi information: tx
load("results/txi/txi_tx.rda")
txi$abundance <- txi$abundance[, !(colnames(txi$abundance) %in% outliers)]
txi$counts <- txi$counts[, !(colnames(txi$counts) %in% outliers)]
txi$length <- txi$length[, !(colnames(txi$length) %in% outliers)]
save(txi, file = "results/txi/txi_tx.rda")
