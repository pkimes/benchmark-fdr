#' Common BenchDesign Design
#' 
#' This function generates a BenchDesign object from the SummarizedBenchmark
#' package that includes several FDR correction methods to be compared in
#' all data analyses. 
#' 
#' Note that before running buildBench() and adding a data object to the output
#' of this function, you'll need to source in any scripts that define wrapper 
#' functions, such as the `scott_fdrreg_hickswrapper` found in 
#' `simulation-helpers.R`.
#'
#' The design assumes the data.frame/list to be benchmarked will have the
#' following columns:
#' - `pval`
#' - `ind_covariate`
#' - `effect_size`
#' - `SE`
#' - `test_statistic`
#' 
#' @return an object of class `BenchDesign` that is initialized with the 
#' benchmarking methods to compare FDR control. These are currently: unadjusted,
#' Bonferroni, Benjamini-Hochberg, Storey's q-value, IHW (with 10 different 
#' alpha values from 0.01 to 0.10), Stephen's ash, Boca-Leek (with 4 different
#' smoothing df values from 2 to 5), Cai's local FDR, and Scott's FDR regression
#' (with two null settings, empirical and theoretical).
#' 
#' @import IHW ashr qvalue swfdr fdrtool FDRreg magrittr SummarizedBenchmark
#' 
#' @author Patrick Kimes

initializeBenchDesign <- function(){
  ## ###########################################################################
  ## check for necessary packages and load them if they aren't present
  ## ###########################################################################
  library(magrittr)
  library(SummarizedBenchmark)
  library(IHW)
  library(ashr)
  library(qvalue)
  library(swfdr)
  library(fdrtool)
  library(FDRreg) 
  
  ## ###########################################################################
  ## define bechmarking design
  ## ###########################################################################

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

  return(bd)
}
