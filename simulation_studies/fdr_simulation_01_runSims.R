#' Alternative Simulations with Uninformative Covariate
#'
#' @author Patrick Kimes

## ##############################################################################
## install packages
## ##############################################################################

## ## from CRAN
## install.packages("ashr")

## ## from Bioconductor
## BiocInstaller::biocLite("IHW")
## BiocInstaller::biocLite("qvalue")
## BiocInstaller::biocLite("swfdr")

## ## from GitHub
## devtools::install_github('jgscott/FDRreg', subdir="R_pkg/")
## devtools::install_github("nignatiadis/IHWpaper")
## devtools::install_github("areyesq89/SummarizedBenchmark", ref = "feat/design-flexibility")


## ##############################################################################
## load packages
## ##############################################################################

## general packages
library("genefilter")
library("dplyr")
library("data.table")
library("ggplot2")
library("magrittr")

## benchmarked methods
library("IHW")
library("ashr")
library("qvalue")
library("swfdr")
library("FDRreg")
library("IHWpaper")

## benchmarking framework
library("SummarizedBenchmark")

## parallelization
library("doParallel")

## helper scripts, e.g. wrappers to FDRreg, clfdr
source("R/simulation-helpers.R")

## set location
projdir <- "/Users/pkimes/workspace/git/projects/project-02-benchmark-fdr"


## ##############################################################################
## define baseline (null) simulation settings
## ##############################################################################

## setting 0: (base) null simulation setting
setting_base <- list(m = 20000,         # number of hypothesis tests
                     pi0 = 1,           # proportion of null hypotheses
                     effect_size = 0,   # expected mean diff of non-null tests
                     n_samples = 20,    # total number of samples
                     n_groups = 2,      # number of groups in contrast
                     informative_ind_covariate = FALSE) # self-explainatory


## ##############################################################################
## define bechmarking design
## ##############################################################################

bd <- BenchDesign()
bd %<>% addMethod("unadjusted", function(x) { x },
                  x = pval)                  
bd %<>% addMethod("bonf",
                  p.adjust,
                  p = pval, method = "bonferroni")
bd %<>% addMethod("bh",
                  p.adjust,
                  p = pval, method = "BH")
bd %<>% addMethod("qvalue",
                  qvalue::qvalue,
                  function(x) { x$qvalues },
                  p = pval)
bd %<>% addMethod("ihw",
                  IHW::ihw,
                  IHW::adj_pvalues,
                  pvalues = pval, covariates = ind_covariate, alpha =  0.1)
bd %<>% addMethod("ashs",
                  ashr::ash,
                  ashr::get_svalue,
                  betahat = effect_size, sebetahat = SE)
bd %<>% addMethod("bl",
                  swfdr::lm_pi0,
                  function(x) { x$pi0 * p.adjust(pval, method = "BH") },
                  pValues = pval, X = ind_covariate, smooth.df = 3)
bd %<>% addMethod("lfdr",
                  clfdr_hickswrapper,
                  unadj_p = pval, groups = IHW::groups_by_filter(ind_covariate, 20))
bd %<>% addMethod("scott",
                  scott_fdrreg_hickswrapper,
                  zscores = qnorm(1 - pval / 2) * sign(effect_size),
                  filterstat = ind_covariate, df = 3, lambda = 0.1,
                  nulltype = 'theoretical')


## ##############################################################################
## define parallelized replications
## ##############################################################################

## number of replications
M <- 10

## set up parallel environment
nCores <- 2
registerDoParallel(cores = nCores)


## ##############################################################################
## Define Simulation Setting 1
##   - no informative covariate
##   - varying null proportion ($\pi_0$)
## ##############################################################################

## create settings w/ varying pi0 (null proportions)
settings <- lapply(seq(.1, 1, by=.1),
                       function(p) { replace(setting_base,
                                             c("pi0", "effect_size"),
                                             c(p, 1.5))
                       })
names(settings) <- paste0("altp0_", 1:10)

## create data.frame of parameters
ptab <- do.call(rbind.data.frame, settings)
ptab$setting <- rownames(ptab)
setnames(ptab, "effect_size", "true_effsize") 

## ##############################################################################
## run simulations and save simulations

## loop over parameter settings
for (idx in seq(settings)) {
    iset <- settings[[idx]]
    
    sblist <- foreach(i = 1:M, .verbose = T) %dopar% {
        ## simulate uninformative indep covariate
        sim_df <- do.call(du_ttest_sim, iset)
        names(sim_df)[which(names(sim_df) == "H")] <- "qvalue"
        
        ## create and process SummarizedBenchmark
        sb <- buildBench(bd, sim_df, truthCol = "qvalue", ptabular = TRUE)        
        sb <- addDefaultMetrics(sb)
        sb <- estimatePerformanceMetrics(sb, alpha=c(0.01, 0.05, 0.1),
                                         addColData = TRUE)
    }
    saveRDS(sblist, file = paste0("data/results-", names(settings)[idx], "-M", M, ".RDS"))
}


## ##############################################################################
## Define Simulation Setting 2
##   - no informative covariate
##   - varying effect size
## ##############################################################################


## ##############################################################################
## run simulations and save simulations



