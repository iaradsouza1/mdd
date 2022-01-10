library(edgeR)
library(ppcseq)
library(tidyverse)
library(magrittr)


# Differentially expressed genes from edgeR -------------------------------

load("results/diff_exp/diff_tx_corrected.rda")

# Sample annotation  ------------------------------------------------------

load("results/important_variables/ann.rda")

# Counts for each gene normalized in TMM ----------------------------------

load("results/txi/txi_tx.rda")

counts <- txi$counts
rownames(counts) <- gsub("\\.\\d+", "", rownames(counts))

##########################################################################

comp <- paste(rep(unique(ann$region), each = 2), unique(ann$gender), sep = "_")

map(comp, function(c) {
    
    cat("comparison: ", c, "\n")
    
    # Diff genes for this comparison - step 1
    diff_tx <- df_res_padj_tx %>% 
        filter(group == c, transcript <= 0.05) %>% 
        unique() %>% 
        pull(txID)
    cat("Step 1 done \n")
    
    # Select information for each region/sex - Step 2 
    ann2 <- ann %>% 
        # tibble::rownames_to_column("sample") %>% 
        mutate(group = paste(region, gender, sep = "_")) %>% 
        filter(group == c) %>% 
        dplyr::select(run, group, ph, rin, phenotype) %>% 
        mutate(phenotype = as.factor(phenotype))
    cat("Step 2 done \n")
    
    # Select counts of diff expressed genes from sex/region comparison - Step 3
    counts2 <- counts %>% 
        as.data.frame() %>% 
        filter(rownames(.) %in% diff_tx)
    cat("Step 3 done \n")
    
    # Transform count data - Step 4
    counts2 <- counts2 %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("tx") %>% 
        gather(key = "run", value = "count", -tx)
    cat("Step 4 done \n")
    
    # Join data and annotation info - Step 5
    counts2 <- inner_join(counts2, ann2, by = "run")
    cat("Step 5 done \n")
    
    # Results from edgeR analysis - Step 6
    df_res <- df_res_padj_tx %>% 
        filter(group == c) %>% 
        dplyr::select(tx = txID, pval = transcript) #Pvalue
    
    ppc_table <- inner_join(counts2, df_res, by = "tx") %>%
        dplyr::select(run, tx, count, phenotype, ph, rin, pval) %>% 
        filter(!is.na(pval))
    cat("Step 6 done \n")
    
    # PPCSEQ - step 7
    counts.ppc <- ppc_table %>%
        mutate(is_significant = pval < 0.05,
               count = as.integer(count+1)) %>%
        identify_outliers(
            formula = ~ phenotype,
            .sample = run, 
            .transcript = tx,
            .abundance = count,
            .significance = pval,
            .do_check = is_significant,
            percent_false_positive_genes = 5, 
            approximate_posterior_inference = FALSE,
            cores = 8
        )
    cat("Step 7 done \n")
    
    return(counts.ppc)
    
}) -> ppc_outliers_tx

names(ppc_outliers_tx) <- comp

# Select only the genes with outliers for further removal
imap_dfr(ppc_outliers_tx, function(x, y) {
    
    idx <- which(x[["tot_deleterious_outliers"]] != 0)
    
    if(length(idx) > 0) {
        
        res <- x %>%
            dplyr::slice(idx) %>%
            dplyr::select(tx, sample_wise_data) %>%
            unnest(cols = sample_wise_data) %>%
            dplyr::select(tx, run, deleterious_outliers) %>%
            filter(deleterious_outliers) %>%
            mutate(group = y) %>%
            dplyr::select(tx, run, group)
        
    }
    
}) -> outliers_samples_dte

save(outliers_samples, file = "results/diff_exp/outliers_samples_dge.rda")
