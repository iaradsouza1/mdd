library(dplyr)
library(purrr)
library(ggplot2)
library(rtracklayer)
library(GenomicFeatures)

# 1. Carregando GTF ---------------------------------------
gtf <- "./data/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz"

gtf_data <- import(gtf)

# 2. Lendo DGE/DTE e pegando os biotipos ---------------------------------------
load("./results/diff_exp/diff_df.rda")

dge_genes <- diff_df %>%
  filter(type == "DGE")

# 2.1. Pegando biotipo dos DGE ---------------------------------------
dge_w_biotype <- gtf_data[, c("gene_id", "gene_biotype")] %>%
  as.data.frame() %>%
  filter(gene_id %in% dge_genes$gene) %>%
  dplyr::select(gene_id, gene_biotype) %>%
  right_join(dge_genes, by = c("gene_id" = "gene")) %>%
  distinct() %>% 
  dplyr::select(gene_id, gene_biotype, group)

readr::write_csv(dge_w_biotype, "results/diff_exp/dge_w_biotype.csv")

# 2.2. Pegando biotipo dos DTE ---------------------------------------
load("./results/diff_exp/diff_tx_corrected.rda")

dte_genes <- df_res_padj_tx %>%
  dplyr::select(txID, transcript, group) %>%
  filter(transcript < 0.01)

dte_w_biotype <-
  gtf_data[, c("transcript_id", "transcript_biotype")] %>%
  as.data.frame() %>%
  filter(transcript_id %in% dte_genes$txID) %>%
  dplyr::select(transcript_id, transcript_biotype) %>%
  right_join(dte_genes, by = c("transcript_id" = "txID")) %>%
  distinct() %>% 
  dplyr::select(transcript_id, transcript_biotype, group)

readr::write_csv(dte_w_biotype, "results/diff_exp/dte_w_biotype.csv")

# 2.3. Pegando biotipo dos DTU ---------------------------------------

dtu_w_biotype <- readr::read_csv("results/ISA/dtu_w_biotype.csv")

# 3. Plotando as porcentagens ---------------------------------------
plot_biotype_bar <- function(data, id_col, n_col) {
  
  id_col <- enquo(id_col)
  n_col <- enquo(n_col)
  
  data %>% 
    ggplot(aes(x = reorder(!!id_col, dplyr::desc(!!n_col)), y = !!n_col)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    coord_flip() +
    labs(
      y = "Porcentagem de Genes",
      x = "Biotipo",
    )
  
}

dge_plot <- dge_w_biotype %>% 
  group_by(gene_biotype) %>% 
  summarise(biotype_n = n() / length(unique(dge_w_biotype$gene_id)) * 100) %>% 
  ungroup() %>% 
  plot_biotype_bar(. , id_col = gene_biotype, n_col = biotype_n)

dte_plot <- dte_w_biotype %>%
  group_by(transcript_biotype) %>%
  summarise(biotype_n = n() / length(unique(dte_w_biotype$transcript_id))* 100) %>%
  ungroup() %>%
  plot_biotype_bar(., id_col = transcript_biotype, n_col = biotype_n)

dtu_plot <- dtu_w_biotype %>%
  group_by(iso_biotype) %>%
  summarise(biotype_n = n() / length(unique(dtu_w_biotype$isoform_id))* 100) %>%
  ungroup() %>%
  plot_biotype_bar(., id_col = iso_biotype, n_col = biotype_n)

ggsave(dge_plot, filename = "results/diff_exp/dge_biotypes.pdf")
ggsave(dte_plot, filename = "results/diff_exp/dte_biotypes.pdf")
ggsave(dtu_plot, filename = "results/diff_exp/dtu_biotypes.pdf")
