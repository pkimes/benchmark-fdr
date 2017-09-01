# Real Data Sets

## Summary of Previous Real Data Analyses

### IHW Real Data Analysis

- ***Application 1: Differential Gene Expression***
  - Raw Data: bulk RNAseq from different mouse strains
  - Input to IHW: 
    - p-values: returned by deseq2
    - covariate: mean of normalized counts for each gene across samples
  - Evaluation metric: number of discoveries
  - Compared with: BH

- ***Application 2: Differential protein abundance***
  - Raw Data: counts from 2666 proteins from quantitative mass-spectometry of yeast subjected to 2 separate treatments: (1) rapamycin and (2) dimethyl sulfoxide
  - Input to IHW: 
    - p-values: differential abundance evaluated by Welch's t-test
    - covariate: number of peptides quantified across all samples for each protein
  - Evaluation metric: number of discoveries
  - Compared with: BH
  
  - ***Application 3: Association between SNPs and histone modifications***
  - Raw Data: SNPs and histone modifications on Chromosome 21
  - Input to IHW: 
    - p-values: 
    - covariate: genomic distance between SNP and ChIP-seq signal
  - Evaluation metric: number of discoveries
  - Compared with: BH


### ASH Real Data Analysis
*summary to be added - e.g. data types, evaluation metrics, code/data availability*


### Boca-Leek Real Data Analysis
*summary to be added - e.g. data types, evaluation metrics, code/data availability*

