# ENRICHMENT ANALYSIS -----------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(patchwork)

# CREATE FUNCTIONS

# Convert IDS
convert_id <- function(data, col) {
  data %>% 
    pull({{col}}) %>% 
    unique() %>% 
    bitr(fromType = "ENSEMBL",
         toType = "ENTREZID",
         OrgDb = "org.Hs.eg.db") %>% 
    mutate(group = data$group[1]) %>%
    mutate(sex = data$sex[1])
}

# Enrichment with enrichGO function
enrich <- function(entrezgenes, col, id) {
  
  res <- enrichGO(gene          = entrezgenes[[col]],
                  OrgDb         = org.Hs.eg.db,
                  keyType       = id,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE) # Uncomment to change ENSEMBL id to SYMBOL
}

# Save plots
create_plots <- function(enrichment, annotation, directory){
    names(enrichment) <- lapply(annotation, function(x){
     unique(as.vector(x$group))
    }) %>% unlist() %>% unname()
    
    if(!dir.exists(directory)) {
      dir.create(directory, recursive = T)
    }
    
    for(i in 1:12){
      if(is.null(enrichment[i][[1]])){
        message(paste(names(enrichment)[i], "is NULL",  sep= ' '))
        next()
      }
      else if(sum(enrichment[[i]]@result$p.adjust < 0.05) == 0){
        message(paste(names(enrichment)[i], "doesn't have enriched terms",  sep= ' '))
        next()
      } else{
        suppressMessages(dotplot(enrichment[[i]], showCategory=20) + scale_color_viridis_c())
        ggsave(filename = paste0(directory, "enrichment_", names(enrichment)[i], ".pdf"),
               width=12, height=8, plot = last_plot())
      }
    }

    return(enrichment)
  }

# Transform results to Dataframe
transform_to_df <- function(enrichment, padj_cutoff = 0.05){
  list_of_df <- list()
  for(i in 1:12){
    if(is.null(enrichment[i][[1]])){
      print(NULL)
    } else{
      enrichment_df <- as_tibble(enrichment[[i]]@result) %>% filter(p.adjust < padj_cutoff) %>% 
      mutate(group = names(enrichment)[i]) %>% list()
      list_of_df <- append(list_of_df, enrichment_df)
    }
  }
  return(list_of_df)
}

# Save dataframes to csv files
save_df <- function(enrichment, file){
  enrichment %>%
    bind_rows() %>%
    write_csv(file=file)
}


# ENRICHMENT OF DIFFERENTIAL GENES ---------------------------------------
load("results/diff_exp/diff_df.rda")

split_genes <- diff_df %>% group_split(group, .keep = T)

converted <- lapply(split_genes, convert_id, col = "gene")

enriched <- lapply(converted, enrich, col = "ENTREZID", id = "ENTREZID")

enriched_ann <- create_plots(enrichment = enriched, directory = 'results/tx_enrich/', annotation = converted)

enriched_df_diff <- transform_to_df(enriched_ann)

save_df(enriched_df_diff, file = 'results/tx_enrich/go_terms.csv')

save(enriched_df_diff, file = "results/tx_enrich/go_terms_tx_by_group.rda")


# ## INTERSECTION ---------------------------------------------------------------
# load("results/intersects/intersects.rda")
# 
# split_genes <- genes_by_sex %>% 
#   group_split(intersect)
# 
# converted <- lapply(split_genes, convert_id, col = "genes")
# 
# enriched <- lapply(converted, enrich, col = "ENTREZID", id = "ENTREZID")
# 
# enriched_ann <- create_plots(enrichment = enriched, directory = 'results/tx_enrich/intersection_by_sex/', annotation = converted)
# 
# enriched_df_by_sex <- transform_to_df(enriched_ann)
# 
# names(enriched_df_by_sex) <- c('female', 'female:male', 'male')
# 
# save_df(enriched_df_by_sex, file = 'results/tx_enrich/intersection_by_sex/go_terms_by_sex.csv')
# 
# save(enriched_df_by_sex, file = "results/tx_enrich/intersection_by_sex/go_terms_by_sex.rda")

# PLOTS -------------------------------------------------------------------

map(enriched_df_diff, ~ {
  if(nrow(.x) > 0) {
    .x %>% 
      mutate(gr = Count/as.numeric(sapply(strsplit(GeneRatio, "/"), "[", 2))) %>% 
      separate(group, into = c("region", "sex"), sep = "_") %>% 
      mutate(sex = str_to_title(sex)) -> df
    
    ggplot(df, aes(y = Description, x = factor(1), size = Count, color = gr)) +
      geom_point() + 
      scale_color_viridis_c(limits = c(0, 0.4)) +
      scale_size_continuous(limits = c(3, 30), breaks = seq(3, 30, 5)) +
      labs(x = "", y = "", color = "Gene ratio") + 
      theme_bw() + 
      labs(title = paste(unique(df$region), unique(df$sex), sep = " ")) + 
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        strip.text = element_text(size = 12),
        panel.spacing.y = unit(2, "lines")
      )
  }
}) -> ls # # List of plots

# Order plots by the number of terms
ls <- ls[sapply(ls, function(x) !is.null(x))]
or <- map(enriched_df_diff, ~ if(nrow(.x) > 0) nrow(.x))
or <- order(unlist(or[sapply(or, function(x) !is.null(x))]), decreasing = T)
ls <- ls[or]

# Plot grid
combined <- wrap_plots(ls) & theme(legend.position = "right")
combined + plot_layout(guides = "collect", heights = c(1, 0.3))

# Save
ggsave("results/tx_enrich/enrichment_diff_regions.pdf", height = 18 , width = 12)