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

## helper digest function
library("digest")

## set location
projdir <- "/n/irizarryfs01/pkimes/projects/project-02-benchmark-fdr"


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
for (idf in 1:5) {
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
                   zscores = qnorm(1 - pval / 2) * sign(effect_size),
                   filterstat = ind_covariate, df = 3, lambda = 0.1,
                   nulltype = 'theoretical')
## Scott's FDR regression w/ empirical null
bd %<>% addBMethod("scott-empirical",
                   scott_fdrreg_hickswrapper,
                   zscores = qnorm(1 - pval / 2) * sign(effect_size),
                   filterstat = ind_covariate, df = 3, lambda = 0.1,
                   nulltype = 'empirical')


## ##############################################################################
## define parallelized replications
## ##############################################################################

## number of replications
M <- 100

## set up parallel environment
nCores <- 10
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


## ##############################################################################
## run simulations and save simulations

## loop over parameter settings
for (idx in seq(settings)) {
    iset <- settings[[idx]]
    
    sblist <- foreach(i = 1:M, .verbose = T) %dopar% {
        ## will break if i > ~2000 (integer seed value too large)
        sim_seed <- (as.integer(Sys.time()) %% 1e6) + (i * 1e6)
        set.seed(sim_seed)
        
        ## simulate data with seed
        sim_df <- do.call(du_ttest_sim, iset)
        names(sim_df)[which(names(sim_df) == "H")] <- "qvalue"

        ## calc data digest
        sim_digest <- digest::sha1(sim_df)

        ## create SummarizedBenchmark
        sb <- buildBench(bd, sim_df, truthCol = "qvalue", ptabular = TRUE)

        ## add simulation information to metadata
        metadata(sb)$sim_func <- du_ttest_sim
        metadata(sb)$sim_parameters <- iset
        metadata(sb)$sim_seed <- sim_seed
        metadata(sb)$sim_digest <- sim_digest
    }
    saveRDS(sblist, file = paste0("data/results-", names(settings)[idx], "-M", M, ".RDS"))
}



## ##############################################################################
## Define Simulation Setting 2
##   - no informative covariate
##   - varying effect size
## ##############################################################################

