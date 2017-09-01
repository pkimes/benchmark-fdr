# Real Data Sets

## Summary of Previous Real Data Analyses

### IHW Real Data Analysis

- ***Application 1: Differential Gene Expression***
  - Raw Data: bulk RNAseq from different mouse strains (Bottomly et al from Recount)
  - Input to IHW: 
    - p-values: returned by deseq2 (default settings, design~strain)
    - covariate: mean of normalized counts for each gene across samples (13 equally sized bins)
    - IHW settings: nfolds=nfolds, cv =ncplits,cv=5
  - Evaluation metric: number of discoveries
  - Compared with: BH

- ***Application 2: Differential protein abundance***
  - Raw Data: counts from 2666 proteins from quantitative mass-spectometry of yeast subjected to 2 separate treatments: (1) rapamycin and (2) dimethyl sulfoxide 
  - Input to IHW: 
    - p-values: differential abundance evaluated by Welch's t-test (supplementary table from [Dephoure, N. et al](https://www.ncbi.nlm.nih.gov/pubmed/22457332))
    - covariate: number of peptides quantified across all samples for each protein (4 equal sized bins)
    - IHW settings: nfolds=nfolds, cv =ncplits,cv=5, regularization parameter selected from a grid of 20 equidistant values
  - Evaluation metric: number of discoveries
  - Compared with: BH
  
  - ***Application 3: Association between SNPs and histone modifications***
  - Raw Data: SNPs and histone modifications on Chromosome 21 from [Grubert, F. et al](https://www.ncbi.nlm.nih.gov/pubmed/26300125)
  - Input to IHW: 
    - p-values: [Matrix eQTL](https://www.ncbi.nlm.nih.gov/pubmed/22492648)
    - covariate: genomic distance between SNP and ChIP-seq signal (10kb bins upto 300kb, 100kb bins upto 1Mb, 10Mb bins for the rest)
    - IHW settings: nfolds=5, \lambda= \inf, cv =ncplits,cv=5, repeated for nominal levels \alpha \subset [0.05,0.1] (grid of 5 equidistant values)
  - Evaluation metric: number of discoveries
  - Compared with: BH, Independent filtering with thresholds set to 10kb, 200kb, 1Mb


### ASH Real Data Analysis
*summary to be added - e.g. data types, evaluation metrics, code/data availability*


### Boca-Leek Real Data Analysis
*summary to be added - e.g. data types, evaluation metrics, code/data availability*

