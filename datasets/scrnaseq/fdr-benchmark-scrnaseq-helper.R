#### helper functions for the scRNAseq analyses

#' compute_wilcoxde_sc
#'
#' @param g_data_sub # the scrnaseq dataframe with gene expression in tpm
#' @param fname #name of RDS file to be saved
#' @param groups #the conditions to be tested
#'
#' @return #none, saves an RDS with computed p-values for each gene along with runtime
#' @export
#'
#' @examples
compute_wilcoxde_sc <- function(g_data_sub,fname,groups){
  
## code borrowed from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R
tmm <- edgeR::calcNormFactors(g_data_sub)
tpmtmm <- edgeR::cpm(g_data_sub, lib.size = tmm * colSums(g_data_sub))
idx <- 1:nrow(tpmtmm)
names(idx) <- rownames(tpmtmm)
timing <- system.time(wilcox_p <- sapply(idx, function(i) {
  wilcox.test(tpmtmm[i, ] ~ groups)$p.value
}))
saveRDS(wilcox_p,fname)
}


compute_scdd_sc <- function(g_cdata_sub,fname1,fname2,fname3){
  Sizes=MedianNorm(g_cdata_sub)
  ng_cdata_sub = Sizes*g_cdata_sub
  saveRDS(Sizes,fname1)
  condition=c(rep(1,(dim(g_cdata_sub)[2])/2),rep(2,(dim(g_cdata_sub)[2])/2))
  rownames(ng_cdata_sub) <- paste0(rownames(ng_cdata_sub), 1:nrow(ng_cdata_sub), sep="")
  colnames(ng_cdata_sub) <- names(condition) <- paste0("Sample",
                                                       1:ncol(ng_cdata_sub), sep="")
  sce <- SingleCellExperiment(assays=list("normcounts"=ng_cdata_sub),
                              colData=data.frame(condition))
  saveRDS(sce,fname2)
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  param=BiocParallel::MulticoreParam(workers=8)
  scDatEx <- scDD(sce, prior_param=prior_param, testZeroes=FALSE)
  ## Notice: 16073 genes have less than 3 nonzero cells per condition.  Skipping these genes.
  
  ##Notice: 16096 genes have less than 3 nonzero cells per condition.  Skipping these genes.
  ##Clustering observed expression data for each gene
  ##Setting up parallel back-end using 6 cores
  saveRDS(scDatEx,fname3)
  
}


compute_MASTtpmde <- function(g_data_sub,groups,fname) {
  ## Acknowledgment: Code borrowed from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTtpm.R
  timing <- system.time({
    sca <- FromMatrix(exprsArray = log2(g_data_sub + 1), 
                      cData = data.frame(wellKey = names(groups), grp = groups))
    zlmdata <- zlm.SingleCellAssay(~groups, sca)
    mast <- lrTest(zlmdata, "groups")
  })
  
  hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
  # list(timing = timing,
  #      res = mast,
       mastdf = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
  saveRDS(mastdf,fname)
}


compute_covariates_scrnaseq <- function(tpm_data_rds){
  g_data_sub=readRDS(tpm_data_rds)
  g_data_sub_med=apply(g_data_sub,1,median)
  g_data_sub_mean=apply(g_data_sub,1,mean)
  g_data_sub_dr=apply(g_data_sub,1,function(x) (length(which(x==0))/length(x)))
  df=data.frame(genes=rownames(g_data_sub),median_exp=g_data_sub_med,mean_exp=g_data_sub_mean,detection_rate=(1-g_data_sub_dr))
  return(df)
}


plot_cov_diagnosis <- function(cov_data,fname){
  png(paste0(fname,"_median_exp_diagnostic.png"),2000,1000,res=100)
  p=strat_hist(cov_data,pvalue ="pval",covariate="median_exp",numQ=5)
  print(p)
  dev.off()
  png(paste0(fname,"_mean_exp_diagnostic.png"),2000,1000,res=100)
  p=strat_hist(cov_data,pvalue ="pval",covariate="mean_exp",numQ=5)
  print(p)
  dev.off()
  png(paste0(fname,"_det_rate_diagnostic.png"),2000,1000,res=100)
  p=strat_hist(cov_data,pvalue ="pval",covariate="detection_rate",numQ=5)
  print(p)
  dev.off()
  png(paste0(fname,"_median_exp_rankscatter.png"),2000,1000,res=100)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="median_exp")
  print(p)
  dev.off()
  png(paste0(fname,"_mean_exp_rankscatter.png"),2000,1000,res=100)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="mean_exp")
  print(p)
  dev.off()
  png(paste0(fname,"_det_rate_rankscatter.png"),2000,1000,res=100)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="detection_rate")
  print(p)
  dev.off()
}
