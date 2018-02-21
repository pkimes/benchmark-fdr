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

print_plot_cov_diagnosis <- function(cov_data){
  print("median_exp")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="median_exp",numQ=5,maxy=10)
  #print(p)
  
  p=strat_hist(cov_data,pvalue ="pval",covariate="median_exp",numQ=3,maxy=10)
  print(p)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="median_exp")
  print(p)
  
  print("mean_exp")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="mean_exp",numQ=5,maxy=10)
  #print(p)
  p=strat_hist(cov_data,pvalue ="pval",covariate="mean_exp",numQ=3,maxy=10)
  print(p)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="mean_exp")
  print(p)
  
  print("detection_rate")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="detection_rate",numQ=5,maxy=10)
  #print(p)
  p=strat_hist(cov_data,pvalue ="pval",covariate="detection_rate",numQ=3,maxy=10)
  print(p)
  p=rank_scatter(cov_data,pvalue ="pval",covariate="detection_rate")
  print(p)
}

cowplot_cov_diagnosis <- function(cov_data){
  #print("median_exp")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="median_exp",numQ=5,maxy=10)
  #print(p)
  
  p1=strat_hist(cov_data,pvalue ="pval",covariate="median_exp",numQ=3,maxy=10)
  p2=strat_hist(cov_data,pvalue ="pval",covariate="mean_exp",numQ=3,maxy=10)
  p3=strat_hist(cov_data,pvalue ="pval",covariate="detection_rate",numQ=3,maxy=10)
  
  p4=rank_scatter(cov_data,pvalue ="pval",covariate="median_exp")
  p5=rank_scatter(cov_data,pvalue ="pval",covariate="mean_exp")
  p6=rank_scatter(cov_data,pvalue ="pval",covariate="detection_rate")
  mean_p=plot_grid(p4,p5,p6, ncol = 3,nrow=1) +draw_label("median_exp", x = 0.15, y = 0.97,fontface="bold") +draw_label("mean_exp", x = 0.45, y = 0.97,fontface="bold")+draw_label("det_rate", x = 0.85, y = 0.97,fontface="bold")
  
  #print("detection_rate")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="detection_rate",numQ=5,maxy=10)
  #print(p)
  #print("mean_exp")
  #p=strat_hist(cov_data,pvalue ="pval",covariate="mean_exp",numQ=5,maxy=10)
  #print(p)
  
  #detr_p=plot_grid(p3, p4, labels = c("mean_exp", "mean_exp"), ncol = 1,nrow=2)
  p<-plot_grid(p1, p2, p3, mean_p, ncol = 1,nrow=4) +draw_label("median_exp", x = 0.5, y = 0.99,fontface="bold") +draw_label("mean_exp", x = 0.5, y = 0.75,fontface="bold")+draw_label("det_rate", x = 0.5, y = 0.5,fontface="bold")
  return(p)
  
  #p<-plot_grid(p1, p2, p3, mean_p, labels = c("median_exp", "mean_exp", "det_rate", "rank_scatter"), ncol = 1,nrow=4),label_x = c(1,1,1,1), label_y = c(0,0,0,0)   
}

plot_covariate_diagnostics <- function(fname,method,h_cov){
  if (method=="wilcox"){
  pval_data=readRDS(fname)
  pval_data_cov=data.frame(genes=names(pval_data),pval=pval_data)
  pval_data_cov=pval_data_cov%>%left_join(h_cov,by=c("genes"))
  p <- cowplot_cov_diagnosis(pval_data_cov)}
  else if (method=="mast"){
    pval_data=readRDS(fname)
    pval_data_cov=data.frame(genes=rownames(pval_data),pval=pval_data$pval)
    pval_data_cov=pval_data_cov%>%left_join(h_cov,by=c("genes"))
    p<-cowplot_cov_diagnosis(pval_data_cov)
  }
  else if (method=="scdd"){
    pval_data=readRDS(fname)
    pval_data_cov=data.frame(genes=rownames(metadata(pval_data)$Genes),pval=metadata(pval_data)$Genes$nonzero.pvalue)
    pval_data_cov$genes=h_cov$genes #scdd changes gene names by appending numeric indices 
    pval_data_cov=pval_data_cov%>%subset(!is.na(pval))
    pval_data_cov=pval_data_cov%>%left_join(h_cov,by=c("genes"))
    p<-cowplot_cov_diagnosis(pval_data_cov)
  }
  return(list(plot=p,data=pval_data_cov))
}



