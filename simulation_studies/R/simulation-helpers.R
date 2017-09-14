#' Simulate Collection of T-Tests
#'
#' Function for simulating a large set of (independent) t-tests. Simulation is performed
#' at the "raw data" level. Test results are returned as a data.frame with one row per
#' test. 
#' 
#' @param m number of hypothesis tests 
#' @param pi0 proportion of true null hypotheses
#' @param effect_size expected mean difference
#' @param n_samples number of total samples, should be multiple of `n_groups` (default = 20)
#' @param n_groups number of groups in comparison (default = 2)
#' @param seed integer seed for random number generator, ignored if NULL (default = NULL) 
#'
#' @return
#' data.frame of test results for `m` simulated data sets, with columns:
#' * `H`: 0/1 indicator whether data simulated under null (0) or alternative (1)
#' * `test_statistic`: t-test statistic
#' * `effect_size`: group mean difference
#' * `pval`: t-test p-value
#' * `ind_covariate`: in this case, it's standard deviation (sd) across all samples 
#' * `SE`: standard error across all samples (sd / sqrt(n))
#'
#' @details
#' This function builds on `IHWpaper::du_ttest_sim` originally written by Nikos Ignatiadis
#' and released under the Artistic-2.0 license. The original function is available at:
#' https://github.com/nignatiadis/IHWpaper
#'
#' @md
#' @importFrom genefilter rowtests
#' @importFrom genefilter rowSds
#' @author Stephanie Hicks
du_ttest_sim <- function(m, pi0, effect_size, n_samples = 20, n_groups = 2, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    m0 <- ceiling(m*pi0)
    false_nulls <- sample(1:m, m - m0)
    
    z_table <- matrix(rnorm(n_samples*m), ncol=n_samples)
    z_table[false_nulls, (n_samples / 2 + 1):n_samples] <- matrix(rnorm(n_samples/2*(m-m0), effect_size), ncol=n_samples/2)

    H <- rep(0, m)
    H[false_nulls] <- 1
    gF <- factor(rep(1:n_groups, each=n_samples/n_groups))
    t_test <- genefilter::rowttests(z_table, gF)
    sds <- genefilter::rowSds(z_table) # pooled var reduces to unpooled var b/c same sample sizes across groups 
    SE <- sds / sqrt(n_samples)
    simDf <- data.frame(H = H, test_statistic = t_test$statistic, effect_size = t_test$dm, 
                        pval = t_test$p.value, ind_covariate = sds, SE = SE)
    return(simDf)
}


#' Calculate Performance Metrics for Adjusted P-Values
#'
#' Function for evaluating the performance of adjusted p-values calculated
#' for a set of t-tests simulated using `du_ttest_sim`. Metrics 
#' calculated for the adjusted p-values are returned in a data.frame with
#' a single row. Each column contains a separate performance metric.
#' 
#' @param sim data.frame output from `du_ttest_sim` containing details of 
#'        simulated hypothesis tests
#' @param adj_p vector of adjusted p-values with length equal to the number
#'        of rows in `sim`
#' @param alpha significance threshold for adjusted p-values
#'
#' @return
#' data.frame with single row summarizing performance of adjusted p-values
#' using various metrics, each in a separate column:
#' * `n_rejects`: total number of rejections
#' * `FWER`: family-wise error rate
#' * `FPR`: false positive rate
#' * `FDP`: false discovery proportion
#' * `power`: power (only returned if at least one H = 1)
#' * `alpha`: alpha input parameter
#'
#' @md
#' @author Stephanie Hicks
calculate_test_stats <- function(sim, adj_p, alpha) {
    ## check that dimensions agree
    stopifnot(length(adj_p) == nrow(sim))
    
    rejected <- adj_p <= alpha
    rjs <- sum(rejected) # number of rejected nulls
    false_rjs <- sum(sim$H == 0 & rejected) # number of false rejected nulls
    rj_ratio <- rjs/nrow(sim)
    FDP <- ifelse(rjs == 0, 0, false_rjs / rjs) # False Discovery Proportion
    power <- ifelse(sum(sim$H) == 0, NA, 
                    sum(sim$H == 1 & rejected) / sum(sim$H == 1)) # power
    FPR <- sum(sim$H == 0 & rejected) / sum(sim$H == 0)
    FWER <- as.numeric(false_rjs > 0)
    if(any(sim$H == 1)){
        df <- data.frame(n_rejects = rjs, FWER=FWER, FPR=FPR, 
                     FDP=FDP, power=power, alpha = alpha) 
    } else{
      df <- data.frame(n_rejects = rjs, FWER=FWER, FPR=FPR, 
                       FDP=FDP, alpha = alpha)
    }
    return(df)
}



#' Simulation Evaluation
#'
#' Wrapper function for computing adjusted p-values (or comparable values) for
#' all methods in benchmkaring study and computing evaluation metrics against
#' expected values at various FDR thresholds.
#' 
#' @param sim data.frame output from `du_ttest_sim` containing details of 
#'        simulated hypothesis tests
#' @param alphas vector of FDR cutoffs between 0, 1
#'
#' @return
#' data.frame with evaluation metrics for each method at the specified alpha
#' thresholds. Should have number of rows equal to `length(alphas) * (n_methods)`,
#' where `n_methods` is the number of different methods (or methods and parameters)
#' in the benchmark study.
#'
#' @md
#' @author Stephanie Hicks
sim_runner <- function(sim, alphas) { 

    ## Bonferroni 
    adj_p <- p.adjust(sim$pval, method="bonferroni")
    df.bonf <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })

    ## BH
    adj_p <- p.adjust(sim$pval, method="BH")
    df.bh <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## qvalue (Storey)
    adj_p <- qvalue::qvalue(p=sim$pval)$qvalues
    df.qvalue <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## IHW (Wolfgang)
    ihw_fdr <- IHW::ihw(pvalues = sim$pval, covariates = sim$ind_covariate, alpha =  0.1)
    adj_p <- IHW::adj_pvalues(ihw_fdr)
    df.ihw <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })

    ## ash (Stephens)
    beta.ash = ashr::ash(betahat = sim$effect_size, sebetahat = sim$SE)
    adj_p <- ashr::get_svalue(beta.ash)
    df.ash <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## swfdr::lm_pi0 (Boca-Leek) **think about multiple covariates** or **correlation structure**
    pi0x <- swfdr::lm_pi0(pValues = sim$pval, X = sim$ind_covariate, smooth.df =3)
    adj_p <- p.adjust(sim$pval, method="BH")*pi0x$pi0
    df.BL <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## locfdr (Cai)
    groupID <- IHW::groups_by_filter(sim$ind_covariate, 20) # stratifies tests based increasing value of an independent covariate
    adj_p <- clfdr_hickswrapper(unadj_p = sim$pval, groups = groupID)
    df.locfdr <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## Scott et al. (2015) (FDR regression (FDRreg) available via GitHub for version 2.0)
    adj_p <- scott_fdrreg_hickswrapper(unadj_p = sim$pval, filterstat = sim$ind_covariate, 
                                       df=3, lambda=0.1, nulltype='theoretical')
    df.scott <- plyr::ldply(alphas, function(a){
        calculate_test_stats(sim=sim, adj_p = adj_p, alpha = a) })
    
    ## Summary table of FDP for each method at each alpha
    return(rbind(data.frame(df.bonf, "method" = "bonferroni"), 
                 data.frame(df.bh, "method" = "bh"), 
                 data.frame(df.qvalue, "method" = "qvalue"), 
                 data.frame(df.ihw, "method" = "ihw"), 
                 data.frame(df.ash, "method" = "ash"), 
                 data.frame(df.BL, "method" = "boca-leek"),
                 data.frame(df.locfdr, "method" = "locfdr"), 
                 data.frame(df.scott, "method" = "scott")))
}


## adapted from IHWpaper::scott_fdrreg()
scott_fdrreg_hickswrapper <- function(unadj_p, filterstat, df=3, lambda=0.01, nulltype = 'theoretical') {
    if (! as.character(packageVersion("FDRreg")) %in% c('0.2.1', '0.2')){
        stop(paste("Benchmarks were run against version 0.2 of FDRreg",
                   "available on github via:",
                   "devtools::install_github(repo= 'jgscott/FDRreg', subdir='R_pkg/')"
                   ))
    }
    
    ## no automated way to choose function space over which we optimize
    ## so we just use bs(df=3) as done in their analysis
    b <- splines::bs(filterstat, df=df)
    
    Xs <- model.matrix( ~  b - 1)
    fdrreg_res <- FDRreg::FDRreg(z=qnorm(unadj_p), features=Xs, nulltype = nulltype,
                                 control=list(lambda = lambda)) # assumption of test statistic follow a standard normal
    adj_p <- fdrreg_res$FDR
    return(adj_p)
}


## adapted from IHWpaper::clfdr()
clfdr_hickswrapper <- function(unadj_p, groups, lfdr_estimation="fdrtool") {
  
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

