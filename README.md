# Repository for scripts used in *Gene- and transcript-level analyses reveal sex-specific trancriptional alterations in prefrontal cortex in Major Depressive Disorder*. 

## Requisites
    
We use Miniconda (or Anaconda) for environment management (See `environment.yml` file). To create the environment:
`conda env create -f environment.yml`

We use R version 4.1.2 and Bioconductor version 3.14.

Here, we use multiple R packages. We use `renv` to track R package versions (see `renv.lock` file).

Clone this repository to your local machine and open the project through `mdd.Rproj` to install required packages. 

## Overall project structure

This project is organized in the following directories:

 - `data/`: Holds the RNA-seq processed data, study metadata, and genome references;  
    - `genome/`: auxiliary files for kallisto quantification.  
        - `Homo_sapiens.GRCh38.97.gtf.gz`: GTF file from human gencode/Ensembl version 97.  
        - `Homo_sapiens.GRCh38.cdna.all.fa.gz`: Human transcriptome gencode/Ensembl version 97.  
    - `kallisto/`: output from kallisto quantification mode (`abundance.tsv`). Quantification files are separated by brain region.  
        - `aINS/`:  
        - `Cg25/`:  
        - `dlPFC/`:  
        - `OFC/`:  
        - `Nac/`:  
        - `Sub/`:  
    - `meta/`:  metadata for BioProject SRP115956.   
        - `SraRunTable.txt` : study metadata  

 - `scripts/`: Holds all the scripts used in the analyses;  
 - `results/`: Holds results from each step analysis.  
  


    







