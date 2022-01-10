library(dplyr)
library(readr)
library(gwasrapidd)

load("results/diff_exp/diff_df.rda")

get_study_variants <- function(study_id) {
  
  variants <- get_variants(study_id = study_id)
  
  var_ensg_ids <- variants@ensembl_ids
  var_variantinfo <- variants@variants %>% 
    select(variant_id, functional_class)
  
  var_ensg_ids %>% 
    left_join(var_variantinfo, by = "variant_id") %>% 
    mutate(study = study_id)
  
}

wray2018 <- 'GCST005839'
howard2018 <- 'GCST005902'
howard2019 <- 'GCST007342'
hyde2016 <- 'GCST006041'

studies <- c(wray2018, howard2018, howard2019, hyde2016)

gwas_mdd <- lapply(studies, get_study_variants) %>% 
  bind_rows()

intersection <- diff_df %>% 
  inner_join(gwas_mdd, by = c("gene" = "ensembl_id")) %>% 
  unique()

intersection %>% 
  write_csv("results/tables/gwas_intersection.csv")

save(intersection, file = "results/diff_exp/gwas_intersections.rda")
