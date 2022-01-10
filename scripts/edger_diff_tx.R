# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS

# GENE

library(edgeR)
library(tidyverse)

# Load metadata -----------------------------------------------------------

load("results/important_variables/ann.rda")

# Set comparisons
comp <- paste(rep(unique(ann$region), each = 2), unique(ann$gender), sep = "_")

# Load gene counts --------------------------------------------------------
load("results/txi/txi_tx.rda")

# Differential expression -------------------------------------------------

# Creating dds object with design matrix including the selected covariate ('rin' and 'ph')
y <- DGEList(counts = txi$counts,
             group = ann$group)

# Design matrix
design <- model.matrix(~ 0 + ph + rin + group, data = ann)
colnames(design)[3:ncol(design)] <- gsub("group", "", colnames(design)[3:ncol(design)])

ct <- makeContrasts(
  aINS_female = aINS_MDD_female-aINS_CTRL_female,
  aINS_male = aINS_MDD_male-aINS_CTRL_male,
  Cg25_female = Cg25_MDD_female-Cg25_CTRL_female,
  Cg25_male = Cg25_MDD_male-Cg25_CTRL_male,
  dlPFC_female = dlPFC_MDD_female-dlPFC_CTRL_female,
  dlPFC_male = dlPFC_MDD_male-dlPFC_CTRL_male,
  Nac_female = Nac_MDD_female-Nac_CTRL_female,
  Nac_male = Nac_MDD_male-Nac_CTRL_male,
  OFC_female = OFC_MDD_female-OFC_CTRL_female,
  OFC_male = OFC_MDD_male-OFC_CTRL_male,
  Sub_female = Sub_MDD_female-Sub_CTRL_female,
  Sub_male = Sub_MDD_male-Sub_CTRL_male, 
  levels = design
)

# Filter low expression genes
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y)

# Estimate dispersions with  
y <- estimateGLMRobustDisp(y, design = design, verbose = T)

# Test
fit <- glmFit(y, design = design)

# Extract results for each comparison

# Summary dataframe
map_df(comp, function(c) {
  
  cat("Comparison: ", c)
  
  lrt <- glmLRT(fit, contrast = ct[, c])
  df <- topTags(lrt, n = Inf, adjust.method = "BH")$table
  if(!is.null(df)) {
    df %>%
      mutate(group = c,
             gene = rownames(df)) %>% 
      dplyr::relocate(gene)
  }
}) -> df_edger_ph_rin_group_tx

# LRT DGElist objects for each comparison
map(comp, function(c) {
  lrt <- glmLRT(fit, contrast = ct[, c])
}) -> lrt_comp
names(lrt_comp) <- comp

save(df_edger_ph_rin_group_tx, lrt_comp, file = "results/diff_exp/edger_tx_rin_ph_diff.rda")


