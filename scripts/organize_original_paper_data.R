
library(readxl)
library(tidyverse)
library(magrittr)

# organize data from original paper 

tm <- read_excel("data/meta/41591_2017_BFnm4386_MOESM3_ESM.xlsx") %>% 
  mutate(sex = "male")

tf <- read_excel("data/meta/41591_2017_BFnm4386_MOESM4_ESM.xlsx") %>% 
  mutate(sex = "female")

df_paper <- tm %>% bind_rows(tf) %>% 
  rename(region = `Brain Region`, gene_name = `Gene name`,
         pvalue = `p-value` ) %>% 
  mutate(
    updown = case_when(
    logFC < 0 ~ "down",
    logFC > 0 ~ "up"
  ),
  region = case_when(
    region == "Anterior Insula" ~ "aINS",
    region == "BA11" ~ "OFC",
    region == "BA25" ~ "Cg25",
    region == "BA8/9" ~ "dlPFC",
    region == "Subic" ~ "Sub",
  ),
  group = paste(region, sex, sep = "_")
    ) %>% 
  select(-`Gene Name`) %>% 
  group_by(group) %>% 
  mutate(padj = p.adjust(pvalue, method = "BH"))

save(df_paper , file = "results/diff_exp/df_original_paper.rda")

# # Load edgeR DE results
# load("~/edger_gene_rin_ph_diff.rda")
# 
# df_edger_ph_rin_group_gene %<>% 
#   mutate(updown = case_when(
#     logFC < 0 ~ "down",
#     logFC > 0 ~ "up"
#   ))
# 
# a <- inner_join(df_edger_ph_rin_group_gene,df_paper, by = c("gene" = "ENSG", "group"))
