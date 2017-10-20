#' Common BenchDesign Design
#' 
#' This script generates a BenchDesign object from the SummarizedBenchmark
#' package that includes several FDR correction methods to be compared in
#' all data analyses. 
#'
#' This script does not load any of the packages need to run these analyses
#' (e.g. the `IHW` package). Packages must be loaded separately in the
#' script where `buildBench` is eventually called.
#'
#' The design assumes the data.frame/list to be benchmarked will have the
#' following columns:
#' - `pval`
#' - `ind_covariate`
#' - `effect_size`
#' - `SE`
#' - `test_statistic`
#' 
#' @author Patrick Kimes

## ##############################################################################
## load packages
## ##############################################################################

library("magrittr")
library("SummarizedBenchmark")

## ##############################################################################
## define bechmarking design
## ##############################################################################

bd <- BenchDesign()
## unadjusted p-values
bd %<>% addBMethod("unadjusted",
                   function(p) { p },
                   p = pval)
## Bonferonni correction
bd %<>% addBMethod("bonf",
                   p.adjust,
                   p = pval, method = "bonferroni")
## Benjamini-Hochberg correction
bd %<>% addBMethod("bh",
                   p.adjust,
                   p = pval, method = "BH")
## Storey's q-value
bd %<>% addBMethod("qvalue",
                   qvalue::qvalue,
                   function(x) { x$qvalues },
                   p = pval)
## IHW (w/ varying alpha threshold) 
for (ia in seq(0.01, 0.10, by=0.01)) {
    bd %<>% addBMethod(paste0("ihw-a", sprintf("%02g", ia*100)),
                       IHW::ihw,
                       IHW::adj_pvalues,
                       pvalues = pval, covariates = ind_covariate,
                       alpha =  UQ(ia))
}
## Stephen's ASH
bd %<>% addBMethod("ashs",
                   ashr::ash,
                   ashr::get_svalue,
                   betahat = effect_size, sebetahat = SE)
## Boca-Leek (w/ varying smoothing degrees of freedom)
for (idf in 2:5) {
    bd %<>% addBMethod(paste0("bl-df", sprintf("%02g", idf)), 
                       swfdr::lm_pi0,
                       function(x) { x$pi0 * p.adjust(pval, method = "BH") },
                       pValues = pval, X = ind_covariate,
                       smooth.df = UQ(idf))
}
## Cai's local FDR
bd %<>% addBMethod("lfdr",
                   clfdr_hickswrapper,
                   unadj_p = pval, groups = IHW::groups_by_filter(ind_covariate, 20))
## Scott's FDR regression w/ theoretical null
bd %<>% addBMethod("scott-theoretical",
                   scott_fdrreg_hickswrapper,
                   zscores = qnorm(1 - pval / 2) * sign(test_statistic),
                   filterstat = ind_covariate, df = 3, lambda = 0.01,
                   nulltype = 'theoretical')
## Scott's FDR regression w/ empirical null
bd %<>% addBMethod("scott-empirical",
                   scott_fdrreg_hickswrapper,
                   zscores = qnorm(1 - pval / 2) * sign(test_statistic),
                   filterstat = ind_covariate, df = 3, lambda = 0.01,
                   nulltype = 'empirical')


