#' Dataset Evaluation (deprecated)
#'
#' This function is no longer used. Please use the SummarizedBenchmark framework 
#' instead! Wrapper function for computing adjusted p-values (or comparable values) and
#' number of rejections at various alphas for
#' all methods in benchmkaring study.
#' 
#' @param dat data.frame containing the following columns: `pval`= p-value for each
#'   test, `ind_covariate` = an independent covariate under the null, 
#'   `effect_size` =  the effect size for each test, `SE` = the standard error 
#'   for each effect. (STILL TO DO - ALLOW MULTIPLE COVARIATES)
#' @param alphas vector of FDR cutoffs between 0, 1
#' @param pvals logical whether to return table of 'adjusted p-values' along with
#'        summary metrics (default = FALSE)
#' @param verbose logical whether to print out progress messages to stdout 
#' 
#' @return
#' If `pvals = FALSE`, a data.frame with evaluation metrics for each method at the specified alpha
#' thresholds. Should have number of rows equal to `length(alphas) * (n_methods)`,
#' where `n_methods` is the number of different methods (or methods and parameters)
#' in the benchmark study.
#' If `pvals = TRUE`, a list with two elements:
#' * `stats`: data.frame containing evaluation metrics
#' * `pvals`: data.frame of the adjusted p-values calculated by
#'    each method, with number of rows equal to `nrows(dat)`.
#'
#' @md
#' @author Stephanie Hicks
#' @author Keegan Korthauer
run_benchmarks <- function(dat, alphas, pvals = FALSE, verbose = TRUE) { 
  
  message("Warning! This function is deprecated! Please don't use this function.",
          " Please use SummarizedBenchmark instead.")
  
  ## keep adjusted p-values
  adj_pset <- list()
  
  ## Bonferroni 
  if(verbose)
    message("Running Bonferonni.")
  adj_p <- p.adjust(dat$pval, method="bonferroni")
  df.bonf <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$bonf <- adj_p
  
  ## BH
  if(verbose)
    message("Running BH.")
  adj_p <- p.adjust(dat$pval, method="BH")
  df.bh <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$bh <- adj_p
  
  ## qvalue (Storey)
  if(verbose)
    message("Running Storey's q.")
  adj_p <- qvalue::qvalue(p=dat$pval)$qvalues
  df.qvalue <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$qvalue <- adj_p
  
  ## IHW (Wolfgang)
  if(verbose)
    message("Running IHW.")
  ihw_fdr <- IHW::ihw(pvalues = dat$pval, covariates = dat$ind_covariate, alpha =  0.1)
  adj_p <- IHW::adj_pvalues(ihw_fdr)
  df.ihw <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$ihw <- adj_p
  
  ## ash (Stephens)
  if(verbose)
    message("Running ash")
  beta.ash = ashr::ash(betahat = dat$effect_size, sebetahat = dat$SE,
                       method="fdr")
  adj_p <- ashr::get_svalue(beta.ash)
  df.ash <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$ashs <- adj_p
  adj_pset$ashq <- ashr::get_qvalue(beta.ash)
  
  ## swfdr::lm_pi0 (Boca-Leek) **think about multiple covariates** or **correlation structure**
  if(verbose)
    message("Running Boca-Leek.")
  pi0x <- swfdr::lm_pi0(pValues = dat$pval, X = dat$ind_covariate, smooth.df =3)
  adj_p <- p.adjust(dat$pval, method="BH")*pi0x$pi0
  df.BL <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$bl <- adj_p
  
  ## locfdr (Cai)
  if(verbose)
    message("Running locfdr.")
  groupID <- IHW::groups_by_filter(dat$ind_covariate, 20) # stratifies tests based increasing value of an independent covariate
  adj_p <- clfdr_hickswrapper(unadj_p = dat$pval, groups = groupID)
  df.locfdr <- plyr::ldply(alphas, function(a){
    rjs <- sum(adj_p <= a)
    data.frame(n_rejects = rjs, alpha = a) 
  })
  adj_pset$lfdr <- adj_p
  
  ## Scott et al. (2015) (FDR regression (FDRreg) available via GitHub for version 2.0)
  # can only run if z-scores are present.
  if("zscore" %in% colnames(dat)){
    if(verbose)
      message("Running Scott.")
    adj_p <- scott_fdrreg_hickswrapper(zscores = dat$zscore, 
                                       filterstat = dat$ind_covariate, 
                                       df=3, lambda=0.1, nulltype='theoretical')
    df.scott <- plyr::ldply(alphas, function(a){
      rjs <- sum(adj_p <= a)
      data.frame(n_rejects = rjs, alpha = a) 
    })
    adj_pset$scott <- adj_p
  }else{
    if(verbose)
      message("Skipping Scott; no Z-scores present")
  }
  
  ## Summary table of FDP for each method at each alpha
  stat_df <- rbind(data.frame(df.bonf, "method" = "bonferroni"), 
                   data.frame(df.bh, "method" = "bh"), 
                   data.frame(df.qvalue, "method" = "qvalue"), 
                   data.frame(df.ihw, "method" = "ihw"), 
                   data.frame(df.ash, "method" = "ash"), 
                   data.frame(df.BL, "method" = "boca-leek"),
                   data.frame(df.locfdr, "method" = "locfdr"), 
                   data.frame(df.scott, "method" = "scott"))
  
  if (pvals) {
    return(list(stats = stat_df,
                pvals = as.data.frame(adj_pset)))
  } else {
    return(stat_df)
  }
}
