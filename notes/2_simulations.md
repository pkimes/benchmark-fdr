# Simulations Settings

## Summary of Previous Simulations

### IWH Simulations
- **Summary:**
  - Simulation studies only consider either BH method or a method with an independent covariate. Goal of simulation studies is to show most covariate adjusted methods result in too many false discoveries, however IHW (1) increases number of discoveries while controlling FDR at alpha level if all null hypotheses are true ("null simulations"), (2) controls FDR at alpha level and increases power as a function of non-null effect sizes. 
  - No modeling assumptions to consider. Simulations are performed either at the level of (1) $p$-values ("null simulations") or (2) $p$-values from two group $t$-tests (normally distributed data with varying effect sizes) where the independent covariate was the pooled variance. 
- **Settings ("null simulations"):** 
  - **methods**:  `IHW:::bh` (BH), `IHW:::lsl_gbh` (GBH with LSL pi0 estimator), `IHW:::tst_gbh` (GBH with TST pi0 estimator), `IHW:::stratified_bh` (stratified BH), `IHW:::clfdr` (wrapper for Cai's local FDR or `locfdr`), `IHWpaper:::scott_fdrreg` (wrapper for `FDRreg`), `IHW:::ddhf` (greedy independent filtering), `IHW::ihw_5fold` (IHW naive with 5 folds)
    - other methods in IHW: `IHW:::storey_qvalue` (`qvalue`), `IHW:::bonf` (bonferonni), `IHW:::gbh` (grouped BH)
    - **settings**: $P_i \sim U[0,1]$, $X_i \sim U[0,1]$, $H_i = 0$ where $i \in (1, m)$ and $m=20000$ with 4000 Monte Carlo replications. FDR control at alpha between 0.01 and 1 (10 equidistant values).  
- **Settings ("varying effect sizes"):** 
  - **methods**:  `IHW:::bh` (BH), `IHW:::lsl_gbh` (GBH with LSL pi0 estimator), `IHW:::tst_gbh` (GBH with TST pi0 estimator), `IHW:::stratified_bh` (stratified BH), `IHW:::clfdr` (wrapper for Cai's local FDR or `locfdr`), `IHWpaper:::scott_fdrreg` (wrapper for `FDRreg`), `IHW:::ddhf` (greedy independent filtering), `IHW::ihw_5fold` (IHW naive with 5 folds)
    - other methods in IHW: `IHW:::storey_qvalue` (`qvalue`), `IHW:::bonf` (bonferonni), `IHW:::gbh` (grouped BH)
  - **settings**: Simulate two groups from normal distribution with means ($\xi_i$) varying from 1 to 2.5 (equidistant grid) and standard deviation of 1. Fraction of true alternatives (pi1) was 0.05. Pooled variance used as the independent covariate. Apply $t$-test and use $p$-values. $m$=20000 and 1000 Monte Carlo replications. 
- **Settings ("size-investing strategy"):** 
  - **idea**: use IHW to shift weight from groups with small $p$-values towards groups with more intermediate $p$-values (since the former will be rejected with a lower weight). Greedy independent filtering, stratified BH, LSL-GBH, TST-GBH, FDRreg cannot apply size investing (and can lose power compared to BH in situations where size investing would be beneficial). 
- **Availability:**
  - Code for simulations, analysis using real data and methods implemented in [IHWpaper are on Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/IHWpaper.html). Use `git clone https://git.bioconductor.org/packages/IHWpaper` to pull down local copy of IHWpaper repository. 


### ASH Simulations
- **Summary:**
  - The simulations focus on two points. **First**, they aim to illustrate how the behavior of the proposed `ash` method and *lfsr* measure, which rely on the "unimodal assumption" (UA), differ from other FDR estimation methods. Naturally, most simulation settings considered follow the UA. **Second**, (and more relevant) they highlight the method's ability to incorporate measurement precision in FDR estimation. Simulations are performed at the level of the test statistics ('beta hats'). 
- **Settings (UA assumption):**
  - **methods**: `ash` (`ash.n`, `ash.u`), `qvalue`, `locfdr` (w/ `nulltype=0`), `mixfdr` (w/ `theonull=TRUE`) (all R packages)
  - **settings**: mixture distributions of test statistics (null + alternative)
    - mixing proportion: pi0 ~ U(0,1)
    - null: point mass at 0
    - alternative: various "g1" densities ("spiky", "near-normal", "flat-top", "skew", "big-normal", "bimodal")
    - estimation error: se=1 for all 
  - **metrics**: accuracy of estimates (*pi0*, *lfdr*, *lfsr*, "g1"), coverage of interval estimates
- **Settings (incorporating precision):**
  - **methods**: `ash`, `qvalue`, `locfdr`
  - **settings**: just one illustrative example; mixture distributions of test statistics
    - mixing proportion: pi0 = 0.5
    - null: point mass at 0
    - alternative: standard normal
    - estimation error: 50% w/ high precision (se=1), 50% w/ low precision (se=10)
- **Availability:**
  - simulations are organized using the author's [dscr](https://github.com/stephens999/dscr/) framework
  - UA simulation code [here](https://github.com/stephenslab/ash/tree/master/code/dsc-shrink)
  - UA simulation parsing [here](https://github.com/stephenslab/ash/blob/master/analysis/index.Rmd#dsc-shrink-the-main-simulation-study-in-the-paper)
  - precision simulation all [here](https://github.com/stephenslab/ash/blob/master/analysis/make_GOODPOOR_figs.Rmd)


### Boca-Leek Simulations
- **Summary:**
  - The goal of the simulations is to demonstrate that their plug-in estimator $\hat{FDR}(x_i)$ controls true FDR while having good power. They compare against other methods that incorporate covariates in FDR estimations (it does not compare against IHW though). 
- **Simulation setting (scenarios)**: They consider 4 scenarios of $\pi_{0}(x_i)$, 
    - (I) A flat function $\pi_{0}(x_i)=0.9$ (no dependence on covariates), 
    - (II) a smooth function in one variable, 
    - (III) a smooth function in one variable within cathegories of a second variable and 
    - (IV) same as in 3, but multiplied by 0.6 to simulate the case where there are less true null hypothesis (or more alternatives). 
- **Simulation setting (Number of tests)**: They consider simulations with $m = 1000$ and $m = 10000$ independent test statistics and added on each 1,000 dependent statistics. For each feature, they sample whether it belongs to the null or alternative distribution based a binomial distribution with probability $\pi_{0}(x_i)$.
- **Simulation setting (Distribution under the alternative)**: Nulls are always assumed to be uniformly distributed between 0 and 1. They considered different distributions for the test statistics under the alternative distribution:
    - $\beta(1, 20)$, 
    - $Norm(\mu, 1)$
    - $T$, 
    - $\chi^{2}$ with 1 df, 
    - $\chi^{2}$ with 4 df, 
    - multivariate normal with either 10 or 20 blocks of correlated tests (with $\rho$'s either .2, .5 or .9) and 
    - T with either 10 or 20 blocks of correlated tests (same $\rho$'s as for the multivariate normal).
- **Simulation**: Given a **scenario**, a **number of tests** and a **distribution for the alternative**, they do 200 runs and estimate the average TPR and FDR. 
- **Methods**: Scott theoretical null, Scott empirical null, Storey and BH.
- **Metrics**: TPR and FDR
- **Availability**: [swdr](https://bioconductor.org/packages/release/bioc/html/swfdr.html) Bioconductor package and several [vignettes](https://github.com/SiminaB/Fdr-regression)
  