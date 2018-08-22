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
#' alpha values from 0.01 to 0.10), Stephens' ash, Boca-Leek (with 4 different
#' smoothing df values from 2 to 5), Cai's local FDR, and Lei and Fithian's
#' AdaPT (with spline GLM).
#' 
#' @import IHW ashr qvalue swfdr fdrtool adaptMT magrittr SummarizedBenchmark
#' 
initializeBenchDesign <- function() {
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
  library(adaptMT)
  library(splines)
  
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
    force(ia)
    bd %<>% addBMethod(paste0("ihw-a", sprintf("%02g", ia*100)),
                       IHW::ihw,
                       IHW::adj_pvalues,
                       pvalues = pval, covariates = ind_covariate,
                       alpha =  UQ(ia))
  }
  ## Stephens' ASH
  bd %<>% addBMethod("ashq",
                     ashr::ash,
                     ashr::get_qvalue,
                     betahat = effect_size, sebetahat = SE)
  ## Boca-Leek (w/ varying smoothing degrees of freedom)
  for (idf in 2:5) {
    force(idf)
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
  ## AdaPT (GLM)
  bd %<>% addBMethod("adapt-glm",
                     adaptMT::adapt_glm,
                     function(x) { x$qvals },
                     pvals = pval,
                     x = data.frame(icov = ind_covariate),
                     pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"),
                     mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"),
                     alphas = 0)
  return(bd)
}


## wrapper to make up for function removed from SummarizedBenchmark package
addDefaultMetrics <- function(sb) {
    sb <- addPerformanceMetric(sb, evalMetric = c("TPR", "FDR", "TNR", "FNR", "rejections"),
                               assay = "qvalue")
    return(sb)
}


## adapted from IHWpaper::clfdr()
clfdr_hickswrapper <- function(unadj_p, groups, lfdr_estimation="fdrtool") {
  
  # Exclude this method if there are fewer than 200 tests within a grouping 
  # (fdrtool is applied separately to each group, and throws a warning in such case)
  if(min(table(groups)) < 200)
    stop("Not enough tests to apply this method. Require at least 200.")
  
  ## estimate local fdr within each stratum first
  lfdr_res <- lfdr_fit(unadj_p, groups, lfdr_estimation=lfdr_estimation)
  lfdrs <- lfdr_res$lfdr
  
  ## now use the rejection rule described in Cai's paper
  
  ## Remark:
  ## When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  ## we get monotonic adjusted p-values as a function of the p-values
  ## This is mainly needed for grenander based lfdrs, with most other
  ## lfdr estimation methods lfdr ties are not a problem usually
  
  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]
  return(adj_p)
}


## helper function for cai
#' @importFrom fdrtool fdrtool
lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){
    
    pvals_list <- split(unadj_p, group)
    
    if (lfdr_estimation == "fdrtool"){
        lfdr_fun <- function(pv) fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
        
    } else if (lfdr_estimation == "locfdr"){
        if (!requireNamespace("locfdr", quietly=TRUE)){
            stop("locfdr package required for this function to work.")
        }
        lfdr_fun <- function(pv) locfdr::locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
    } else {
        stop("This lfdr estimation method is not available.")
    }
    
    lfdr_list <- lapply(pvals_list, lfdr_fun)
    lfdrs <- unsplit(lfdr_list, group)
    
    fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
    fit_obj
}
