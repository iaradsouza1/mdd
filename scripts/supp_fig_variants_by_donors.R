
# Create supplementary figures 9 and 10 -----------------------------------

# This script uses the Metadata.csv file, which summarises the relationship among
# samples from each brain region and the donors. This table was reorganized based
# on the data provided by Labonté et al. (2017) paper. 

library(tidyverse)
library(vcfR)
library(RColorBrewer)
library(ComplexHeatmap)

# Import table with samples ids -------------------------------------------

df <- read_csv("data/meta/Metadata.csv", col_types = paste(rep("c", 20),collapse = ""))
df %>% 
  gather(key = "region", "rin", mvPFC:vSUB) %>% 
  mutate(Gender = tolower(Gender),
         region = case_when(
           region == "mvPFC" ~ "Cg25",
           region == "vSUB" ~ "Sub",
           region == "Nac" ~ "Nac",
           region == "OFC" ~ "OFC",
           region == "aINS" ~ "aINS",
           region == "dlPFC" ~ "dlPFC"
         ),
         pH = as.character(round(as.numeric(pH),2)),
         `Alcohol (y/n)` =  tolower(`Alcohol (y/n)`),
         `Alcohol (y/n)` = ifelse(is.na(`Alcohol (y/n)`), "NA",`Alcohol (y/n)`)) %>% 
  select(Samples, Gender, `Alcohol (y/n)`,`Cause of death`, region, pH, Age, rin, Group) %>% 
  unique() -> df2

# Import data from SRA ----------------------------------------------------

ann <- read.csv("data/meta/SraRunTable.txt", stringsAsFactors = F, header = T)
colnames(ann) <- tolower(colnames(ann))
ann <- ann %>% 
  filter(organism == "Homo sapiens") %>% 
  mutate(
    tissue = case_when(
      tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
      tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
      tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
      tissue == "Anterior Insula (aINS)" ~ "aINS",
      tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
      tissue == "Subiculum (Sub)" ~ "Sub"
    ),
    ph = as.character(round(as.numeric(ph),2)),
    age = as.character(age),
    rin = as.character(signif(rin)),
    alcool = tolower(alcool),
    alcool =  ifelse(is.na(alcool), "NA", alcool)
  ) %>% 
  select(run, gender, alcool, region = tissue, biosample, ph, age, rin, cause_of_death, phenotype) %>% 
  unique() -> ann2

# -------------------------------------------------------------------------

inner_join(ann2, df2, 
           by = c("age" = "Age",
                  "region",
                  "phenotype" = "Group",
                  "ph" = "pH",
                  "cause_of_death" = "Cause of death",
                  "gender" = "Gender",
                  "rin",
                  "alcool" = "Alcohol (y/n)")) -> ann_ind


# Create heatmaps ---------------------------------------------------------

# Plot heatmap of the number of samples by individual
ann_ind %>% 
  group_by(Samples, region) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  spread(key = Samples, value = n) %>% 
  as.data.frame() -> tmp_m

# temporary matrix that holds samples distribution
rownames(tmp_m) <- tmp_m$region
tmp_m$region <- NULL
tmp_m <- t(tmp_m)
tmp_m[is.na(tmp_m)] <- 0
tmp_m <- cbind(tmp_m, rowSums(tmp_m))
colnames(tmp_m)[7] <- "soma"

tmp_m <- tmp_m[order(tmp_m[,"soma"], decreasing = T),]

dir <- "data/vcfs/"

# Plot variants distribution ----------------------------------------------

# Read vcfs 
# VCFs were filtered based on variants that intersected with the ones in Table 1. 
map_dfr(ann_ind$run, function(x) {
  
  df <- read.vcfR(paste0(dir, x, "_intersect_gwas.vcf"))
  print(x)
  if(nrow(df@fix) != 0) {
    res <- as.data.frame(df@fix)
    res$run <- x
    return(res)
  }
}) -> vcf_variants

# Merge variants and individuals info
vcf_variants %>% inner_join(ann_ind %>% select(run, region, Samples), by = "run") %>% 
  mutate(gen = ifelse(sapply(strsplit(vcf_variants$INFO, ";"), "[[", 1) == "AC=1", "Het", "Hom")) -> vcf_variants2

vcf_variants2 %>% 
  select(ID, run, region) %>%
  count(ID, region) %>% 
  spread(region, n) %>% 
  mutate(across(everything(), ~ ifelse(is.na(.), 0, . ))) %>%
  rowwise() %>% 
  mutate(total = sum(c_across(aINS:Sub))) %>% 
  arrange(desc(total)) %>% 
  select(-total) %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix()-> vcf_variants3


# Plot variants by region and by sample -----------------------------------

# Merge info and filter control samples
vcf_variants %>% inner_join(ann_ind %>% select(run, Samples, region, gender, phenotype), by = "run") %>% 
  filter(phenotype == "CTRL") -> vcf_variants2

# Transform dataframe 
vcf_variants2 %>% 
  select(ID, region, run, Samples, region, phenotype) %>%
  unique() %>% 
  spread(key = Samples, value = run) %>% 
  arrange(phenotype, ID, region) %>% 
  mutate(across("128":"74", ~ ifelse(is.na(.), 0, 1))) -> vcf_variants3

# Create matrix to plot and the annotation dataframe
m_ctrl <- vcf_variants3 %>% 
  filter(phenotype == "CTRL") %>% 
  select(-c(1:3))
m_ctrl_ann <- vcf_variants3[,1:3] %>% 
  filter(phenotype == "CTRL")

# Merge info and filter MDD samples
vcf_variants %>% inner_join(ann_ind %>% select(run, Samples, region, gender, phenotype), by = "run") %>% 
  filter(phenotype == "MDD") -> vcf_variants2

vcf_variants2 %>% 
  select(ID, region, run, Samples, region, phenotype) %>%
  unique() %>% 
  spread(key = Samples, value = run) %>% 
  arrange(phenotype, ID, region) %>% 
  mutate(across("105":"93", ~ ifelse(is.na(.), 0, 1))) -> vcf_variants3

m_mdd <- vcf_variants3 %>% 
  filter(phenotype == "MDD") %>% 
  select(-c(1:3))
m_mdd_ann <- vcf_variants3[,1:3] %>% 
  filter(phenotype == "MDD")

# Create heatmaps
pdf("~/Área de trabalho/vcfs/variants_by_sample_by_region.pdf", height = 12, width = 10)
Heatmap(m_ctrl, cluster_rows = F, cluster_columns = F, column_title = "CTRL",
        row_title_rot = 0, col = c("white", "red"),
        rect_gp = gpar(col = "grey60", lwd = 1), row_gap = unit(5, "mm"),
        row_split = m_ctrl_ann$ID, row_labels = m_ctrl_ann$region)

Heatmap(m_mdd, cluster_rows = F, cluster_columns = F, column_title = "MDD",
        col = c("white", "red"), row_gap = unit(5, "mm"),
        rect_gp = gpar(col = "grey60", lwd = 1), row_title_rot = 0,
        row_split = m_mdd_ann$ID, row_labels = m_mdd_ann$region)
dev.off()

# Save table with all results
vcf_variants %>% inner_join(ann_ind %>% select(run, Samples, region, gender, phenotype), by = "run") %>% 
  select(-QUAL, -FILTER, -INFO, -POS) %>% 
  relocate(Samples:phenotype, .before = "CHROM") %>% 
  arrange(Samples, ID) %>% 
  write_csv(file = paste(dir, "variants_gwas_by_individual.csv"), quote = "none")
