# Simulations Settings

## Summary of Previous Simulations

### IWH Simulations
*summary to be added - e.g. settings, evaluation metrics, availability*


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
*summary to be added - e.g. settings, evaluation metrics, availability*

