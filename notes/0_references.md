# References

## Core Methods
- :fire: **IHW (Independent Hypothesis Weighting)** :fire:
  - :books:**Reference:**
    - [Ignatiadis N, Klaus B, Zaugg JB, and Huber W. (2016). "Data-driven hypothesis weighting increases detection power in genome-scale multiple testing." Nature Methods, 13(7):577-580.](https://www.ncbi.nlm.nih.gov/pubmed/27240256)
  - :computer:**Code:**
    - `IHW` R package [(Bioconductor link)](https://bioconductor.org/packages/release/bioc/html/IHW.html)
    - code from analysis [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- :fire: **ASH (Adaptive SHrinkage)** :fire:
  - :books:**Reference:**
    - [Stephens M. (2017). "False discovery rates: a new deal." Biostatistics, 18(2):275-294.](https://www.ncbi.nlm.nih.gov/pubmed/27756721)
  - :computer:**Code:**
    - `ashr` R package [(GitHub link)](https://github.com/stephens999/ashr), [CRAN link](https://cran.r-project.org/web/packages/ashr/index.html)
    - code from analysis [(GitHub link)](https://github.com/stephenslab/ash)
- :fire: **Boca-Leek** :fire:
  - :books:**Reference:**
    - [Boca S, and Leek J. (2017). "A direct approach to estimating false discovery rates conditional on covariates." bioRxiv preprint.](http://www.biorxiv.org/content/early/2017/07/25/035675)
  - :computer:**Code:**
    - `swfdr` R package [(Bioconductor link)](https://bioconductor.org/packages/release/bioc/html/swfdr.html) applies the Boca-Leek regression framework to estimate the Science-Wise False Discovery Rate from [Jager & Leek (2013)](https://arxiv.org/abs/1301.3718)
    - code from analysis of simulated and real data in the Boca-Leek paper: [(GitHub link)](https://github.com/SiminaB/Fdr-regression)

## Relevant Methods
- **Benjamini and Hochberg (BH)**
  - :books:**Reference:**
    - [Benjamini and Hochberg. (1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological), 57(1): 289-300.](http://www.stat.purdue.edu/~doerge/BIOINFORM.D/FALL06/Benjamini%20and%20Y%20FDR.pdf)
  - :computer:**Code:**
    - implemented in the `p.adjust()` function in `stats` R package
- **Weighted BH** (basis of IHW)
  - :books:**Reference:**
    - [Roeder K, Devlin B, Wasserman L. Improving power in genome-wide association studies: weights tip the scale. Genet Epidemiol. 2007 Nov;31(7):741-7.](http://www.stat.cmu.edu/~roeder/publications/rdw2007.pdf)
  - :computer:**Code:**
    - implemented in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **Greedy Independent Filtering** (compared in IHW paper)
  - :books:**Reference:**
    - [Bourgon R, Gentleman R and Huber W. (2010) "Independent filtering increases detection power for high-throughput experiments." Proceedings of the National Academy of Sciences. 107:9546–9551.](http://www.pnas.org/content/107/21/9546.long)
  - :computer:**Code:**
    - implemented in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **LSL-GBH, TST-GBH** (compared in IHW paper)
  - :books:**Reference:**
    - [Hu JX, Zhao H, and Zhou HH. (2010). "False discovery rate control with groups." Journal of the American Statistical Association, 105(491):1215-1227.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3175141/)
  - :computer:**Code:**
    - implemented in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **SBH** (compared in IHW paper)
  - :books:**Reference:**
    - [Sun L, Craiu RV, Paterson AD, and Bull SB. (2006). "Stratified false discovery control for large‐scale hypothesis testing with application to genome‐wide association studies." Genetic Epidemiology, 30(6):519-530.](https://www.ncbi.nlm.nih.gov/pubmed/16800000)
  - :computer:**Code:**
    - implemented in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **Clfdr (Cai's local FDR)** (compared in IHW paper)
  - :books:**Reference:**
    - [Cai TT, and Sun W. (2009). "Simultaneous testing of grouped hypotheses: Finding needles in multiple haystacks." Journal of the American Statistical Association, 104(488):1467-1481.](https://pdfs.semanticscholar.org/b757/c6ca12e6db25a258d1078994e3ad96e18f06.pdf)
  - :computer:**Code:**
    - implemented in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **FDRreg** (compared in IHW paper, swfdr paper)
  - :books:**Reference:**
    - [Scott JG, Kelly RC, Smith MA, Zhou P, and Kass RE. (2015). "False discovery rate regression: an application to neural synchrony detection in primary visual cortex." Journal of the American Statistical Association, 110(510):459-471.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4743052/)
  - :computer:**Code:**
    - `FDRreg` R package [(GitHub link)](https://github.com/jgscott/FDRreg)
    - wrapped in `IHWpaper` [(GitHub link)](https://github.com/nignatiadis/IHWpaper/) [(Bioconductor link)](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html)
- **local false discovery rate** (compared in ASH paper)
  - :books:**Reference:**
    - [Efron B. (2008). "Microarrays, empirical Bayes and the two-groups model." Statistical Science, 23(1):1-22.](http://projecteuclid.org/download/pdfview_1/euclid.ss/1215441276)
  - :computer:**Code:**
    - `locfdr` R package [(CRAN link)](https://cran.r-project.org/web/packages/locfdr/index.html)
- **mixture false discovery rate** (compared in ASH paper)
  - :books:**Reference:**
    - [Muralidharan O. (2010). "An empirical Bayes mixture method for effect size and false discovery rate estimation." The Annals of Applied Statistics, 4(1):422-438.](http://projecteuclid.org/download/pdfview_1/euclid.aoas/1273584461)
  - :computer:**Code:**
    - `mixfdr` R package [(CRAN link, archived)](https://cran.r-project.org/web/packages/mixfdr/index.html)






