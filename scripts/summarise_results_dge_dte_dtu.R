
# Load DGE/DTE/DTU results ------------------------------------------------

library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(biomaRt)
library(magrittr)

# Dictionary of IDs (ensembl_transcript_id, ensembl_gene_id, and gene name) ----
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
tx2gene <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
                 mart = ensembl)

tx2gene$hgnc_symbol <- ifelse(tx2gene$hgnc_symbol == "", NA, tx2gene$hgnc_symbol)
tx2gene$ensembl_gene_id <- ifelse(tx2gene$ensembl_gene_id == "", NA, tx2gene$ensembl_gene_id)
tx2gene$ensembl_transcript_id <- ifelse(tx2gene$ensembl_transcript_id == "", NA, tx2gene$ensembl_transcript_id)

tx2gene <- tx2gene[complete.cases(tx2gene), ]

rm(ensembl)

# DGE and DTE results ----
load("results/diff_exp/df_res_dge_dte.rda")

df_res_padj_gene_out_filtered %<>% 
  inner_join(
    tx2gene %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% unique(), 
    by = c("gene" = "ensembl_gene_id"))

df_res_padj_tx_out_filtered %<>%
  inner_join(
    tx2gene %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>% unique(), 
    by = c("geneID" = "ensembl_gene_id")) %>% 
  unique()

# DTU results ----

# Load 2 step object from ISA
files <- list.files("results/ISA/", pattern = "pass2", recursive = T, full.names = T)
isa_df <- map_dfr(files, ~ {
  load(.x)
  SwitchList_2$isoformFeatures
})

# Conditions to filter results:
condition_1_male <- grepl("CTRL_male", isa_df$condition_1)
condition_2_male <- grepl("MDD_male", isa_df$condition_2)
condition_1_female <- grepl("CTRL_female", isa_df$condition_1)
condition_2_female <- grepl("MDD_female", isa_df$condition_2)

isa_df <- isa_df[(condition_1_male & condition_2_male) | 
                   (condition_1_female & condition_2_female),]

# Organize df
isa_df %<>%
  mutate(gene_id = gsub("\\.\\d+", "", gene_id),
         isoform_id = gsub("\\.\\d+", "", isoform_id)) %>% 
  filter(isoform_switch_q_value <= 0.05) %>% 
  dplyr::select(isoform_id, gene_id, condition_1) %>% 
  inner_join(
    tx2gene %>% dplyr::select(ensembl_transcript_id, ensembl_gene_id) %>% unique(),
    by = c("isoform_id" = "ensembl_transcript_id")
  ) %>% 
  mutate(condition_1 = gsub("_CTRL", "", condition_1)) %>% 
  dplyr::select(gene = ensembl_gene_id, tx = isoform_id, hgnc_symbol = gene_id, group = condition_1) %>% 
  unique()

# Filter out transcripts that were identified as outliers by ppcseq analysis ----

load("results/diff_exp/outliers_samples_dte.rda")
isa_df %<>% anti_join(outliers_samples_dte, by = "tx")

# Combine all 3 analyses --------------------------------------------------

dge <- df_res_padj_gene_out_filtered %>% 
  dplyr::select(gene, hgnc_symbol, group) %>% 
  unique() %>% 
  mutate(type = "DGE")
dte <- df_res_padj_tx_out_filtered %>% 
  dplyr::select(gene = geneID, hgnc_symbol, group) %>% 
  unique() %>% 
  mutate(type = "DTE")
dtu <- isa_df %>% 
  dplyr::select(gene, hgnc_symbol, group) %>% 
  unique() %>% 
  mutate(type = "DTU")

# Main dataframe used in downstream analysis
diff_df <- purrr::reduce(list(dge, dte, dtu), bind_rows)

# Function to count genes in each analysis
count_genes <- function(diff_df, col_filter = "type", col_value, by, id) {
  
  diff_df %>% 
    dplyr::filter(.data[[col_filter]] == col_value) %>% 
    dplyr::select({{id}}, {{by}}) %>% 
    unique() %>%
    group_by(.data[[by]]) %>% 
    summarise(n_genes = n())
  
}

# Number of altered genes in each analysis
count_genes(diff_df, col_filter = "type", col_value = "DGE", by = "group", id = "gene") 
count_genes(diff_df, col_filter = "type", col_value = "DTE", by = "group", id = "gene")
count_genes(diff_df, col_filter = "type", col_value = "DTU", by = "group", id = "gene")

# Number of genes in each analysis
diff_df %>% 
  dplyr::select(type, gene) %>% 
  unique() %>% 
  dplyr::count(type)

save(diff_df, file = "results/diff_exp/diff_df.rda")
  

