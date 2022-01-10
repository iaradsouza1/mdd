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
 - `scripts/`: Holds all the scripts used in the analyses;
 - `results/`: Holds results from each step analysis.

 In each directory, there is detailed description of files. 






