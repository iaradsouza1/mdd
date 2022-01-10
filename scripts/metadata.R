# Based on the results from rank_variables.R analysis, we decided to include only 
# the first 2 covariates: 'rin' and 'ph'.

library(dplyr)

# Read and select the covariates
ann <- read.csv("data/meta/SraRunTable.txt", stringsAsFactors = F, header = T)
colnames(ann) <- tolower(colnames(ann))
ann <- ann %>% 
  dplyr::select(run, ph, rin, phenotype, gender, tissue, organism) %>% 
  filter(organism == "Homo sapiens") %>% 
  mutate(region = case_when(
    tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
    tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
    tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
    tissue == "Anterior Insula (aINS)" ~ "aINS",
    tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
    tissue == "Subiculum (Sub)" ~ "Sub"
  ), 
    group = paste(region, phenotype, gender, sep = "_")) %>% 
  dplyr::rename(sample_id = run) %>% 
  dplyr::select(sample_id, ph, rin, phenotype, gender, region, group)

# Scale ph and rin
ann$ph <- scale(ann$ph)
ann$rin <- scale(ann$rin)
rownames(ann) <- ann$sample_id
ann$sample_id <- NULL

# Save (this metadata will be used in all future models)
save(ann, file = "results/important_variables/ann.rda")

  

