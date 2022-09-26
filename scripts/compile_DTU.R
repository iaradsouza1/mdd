library(dplyr)
library(purrr)
library(IsoformSwitchAnalyzeR)

# 1. Get all results from ISA ------

files <-
  list.files(
    "results/ISA/",
    pattern = "pass2",
    recursive = T,
    full.names = T
  )
isa_df <- map_dfr(files, ~ {
  load(.x)
  SwitchList_2$isoformFeatures
})

# 2. Filter results by conditions and save full result -----

condition_1_male <- grepl("CTRL_male", isa_df$condition_1)
condition_2_male <- grepl("MDD_male", isa_df$condition_2)
condition_1_female <- grepl("CTRL_female", isa_df$condition_1)
condition_2_female <- grepl("MDD_female", isa_df$condition_2)

isa_df <- isa_df[(condition_1_male & condition_2_male) |
                   (condition_1_female & condition_2_female), ]

save(isa_df, file = "results/ISA/DTU_df.rda")

# 3. Get dataframe with biotypes ----

dtu_w_biotypes <- isa_df %>%
  mutate(
    gene_id = gsub("\\.\\d+", "", gene_id),
    isoform_id = gsub("\\.\\d+", "", isoform_id)
  ) %>%
  filter(isoform_switch_q_value <= 0.05) %>%
  mutate(group = gsub("_CTRL", "", condition_1)) %>%
  dplyr::select(isoform_id, iso_biotype, group)

readr::write_csv(dtu_w_biotypes, file = "results/ISA/dtu_w_biotype.csv")
