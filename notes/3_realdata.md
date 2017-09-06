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
- none presented

### Boca-Leek Real Data Analysis
- **Summary:** A meta-analysis GWAS of BMI was used to estimate FDR conditional on sample size and MAF of each SNP. The idea is that increased sample size or MAF closer to 0.5 increases the power to detect associations of SNPs with BMI
- **Dataset:** 
  - GWAS for BMI from the Genetic Investigation of ANthropometric Traits (GIANT) consortium [(Locke et al., 2015)](https://www.ncbi.nlm.nih.gov/pubmed/25673413)
  - Meta-analysis with 339,224 total samples with BMI and 2,555,510 SNPs measured (not every SNP was measured in each sample)
  - Covariates of interest are (1) the sample size for each SNP, and (2) the minor allele frequency (MAF)
- **Input:**
  - p-values (for BL, BH, and Storey), sample size, MAF, and Z-scores (for Scott) from [http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files]( http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)
- **Model Matrix:**
  - SNP-specific sample size included as a natural cubic spline with 5 df (to >allow for sufficient flexibility)
  - MAF categorized into three bins (tertiles)
- **Compared with:** 
  - Including covariates: Scott (2015) with theoretical null (denoted T) and Scott (2015) with empirical null (denoted E)
  - Not including covariates: Storey (2002) and Benjamini-Hochberg (1995)
- **Settings:**
  - FDR level 0.05
  - BL: `smooth.df=3`, `smooth.df=3` (this means lambda is taken as the smoothed value at 0.95 using a cubic spline with 3 df)
  - Scott: `lambda=1`, `nulltype = 'theoretical'` (T) and `nulltype = 'empirical'` (E)
- **Results:**
  - p-values decrease as sample size increases
  - Scott (2015) with theoretical null found substantially more discoveries than BL
  - Scott (2015) with empirical null found substantially fewer discoveries than BL
  - BH and Storey found slightly fewer discoveries than BL
- **Availability:**
  - Code for analysis using real data is available here [(GitHub link)](https://github.com/SiminaB/Fdr-regression)


