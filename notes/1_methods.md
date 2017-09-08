# Methods

## Core Methods
- :fire: **IHW (Independent Hypothesis Weighting)** :fire:
  - :bulb:**Idea:** Method uses a data-driven approach to calculate weights using independent covariates and then applies a group-weighted BH method
  - **Features** [Size investing](https://projecteuclid.org/euclid.aos/1297779856)
  - **Assumptions:**
    - covariate must be independent of p-value under null
    - tests are independent under null 
  - :+1:**Good:** Improves power and controls FDR at alpha level
  - :-1:**Bad:** does not maintain FDR when (1) insufficient power to detect false hypotheses (2) all hypotheses are true
- :fire: **ASH (Adaptive SHrinkage)** :fire:
  - :bulb:**Idea:** 
  - :+1:**Good:**
  - :-1:**Bad:**
- :fire: **Boca-Leek** :fire:
  - :bulb:**Idea:** 
    - Estimates FDR conditional on covariates in a multiple testing framework by (1) estimating the null proportion of hypotheses with logistic regression and (2) multiplying the BH-adjusted p-values by the estimated null proportion.
    - Assumes that the p-values are independent of the covariates conditional on the null or alternative. In other words, the probability that a test statistic or p-value is drawn from the null or non-null distribution depends on the covariates (but not the realized value of the test statistic or p-value itself).
  - :+1:**Good:**
    - FDR control can be maintained even when tests are moderately correlated
    - If hypotheses are independent, can use bootstrapping to obtain CIs for the null proportion of hypotheses
    - Increasing the number of tests can lead to improvement in FDR control
  - :-1:**Bad:**
    - Improvement over Storey (2002) is minor
    - Scott (2015) is superior when test statistics are normally distributed
    - Does not control FDR when hypotheses are highly correlated
    - Requires tuning of $\lambda$ parameter and careful specification of covariates (smoothing splines? categorization?)

## Relevant Methods
- **Benjamini and Hochberg (BH)**
  - :bulb:**Idea:** Uses *p*-values. Assumes independence of tests (or assumes all hypotheses are exchangeable), but has been shown to be robust even with correlated tests. 
  - :+1:**Good:** Easy to use. More powerful than FWER-based methods
  - :-1:**Bad:** Sub-optimal power (and can't prioritize tests) when the individual tests differ in statistical properties such as sample size, true effect size, signal-to-noise ratio or prior probability of being false (see [Ignatiadis et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27240256)) 
    - To increase power while controlling FDR, can use independent covariates to prioritize tests 
- **Weighted BH** (basis of IHW)
  - :bulb:**Idea:** Associate each test with a non-negative weight such that weights averahe to 1. Hypotheses with higher weights are prioritized. 
  - :+1:**Good:** Controls FDR 
  - :-1:**Bad:** Weights must be prespecified, are independent of data
- **Greedy Independent Filtering** (compared in IHW paper)
  - :bulb:**Idea:** Filter all p-values using an independent covariate ($X$) such that $X$ is less than some threshold $x$. Greedy if researcher tests all possible thresholds and picks one that maximizes number of discoveries. However, greedy approaches do not control for FDR at alpha level. 
  - :+1:**Good:** Controls FDR at alpha level IF researcher ONLY uses a pre-specified threshold
  - :-1:**Bad:** Not automated, very subjective, can lead to *p*-hacking
- **LSL-GBH, TST-GBH** (compared in IHW paper)
  - :bulb:**Idea:** 
  - :+1:**Good:**
  - :-1:**Bad:**
- **SBH** (compared in IHW paper)
  - :bulb:**Idea:** Use a categorial (or binned continuous) covariate to stratify tests, apply BH within each strata, combine significant tests
  - :+1:**Good:** 
  - :-1:**Bad:** Loss of FDR control for the null tests because different strata are treated equally
    - To increase power, different strata should be prioritized differently
- **Clfdr (Cai's local FDR)** (compared in IHW paper)
  - :bulb:**Idea:** 
  - :+1:**Good:**
  - :-1:**Bad:**
- **FDRreg** (compared in IHW paper, swfdr paper)
  - :bulb:**Idea:** 
  - :+1:**Good:**
  - :-1:**Bad:**
- **local false discovery rate** (compared in ASH paper)
  - :bulb:**Idea:** local FDR is estimated separately within each group and then estimates are grouped together
  - :+1:**Good:** 
  - :-1:**Bad:** FDR isn't controlled at alpha level if alternative distributions across groups are different
- **mixture false discovery rate** (compared in ASH paper)
  - :bulb:**Idea:** 
  - :+1:**Good:**
  - :-1:**Bad:**






