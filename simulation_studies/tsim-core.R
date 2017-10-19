#' T-Test Simulations for FDR Assessment
#' 
#' This script is should be run from the command line with the following
#' four ordered arguments:
#'
#' 1. M: integer number of simulation replications
#' 2. ncores: integer number of computing cores to use for parallelization
#' 3. setting_vparam: simulation setting (variable parameter)
#' 4. setting_icparam: simulation setting (informative covariate parameter)
#'
#' For a specified simulation settings ('setting_vparam' + 'setting_icparam'),
#' this script applied FDR correction to results from a collection of independent
#' t-tests. For more deetails on the simulation settings, see the `define_settings`
#' function defined in `sim-core-settings.R`. 
#'
#' This script works with the following other scripts:
#' - `R/du_tsim.R`: `du_tsim`
#'     This function is used for actually simulating t-tests.
#' - `sim-core-settings.R`: `define_settings`
#'     This function is used for defining the settings to be passed to
#'     `du_tsim` for the simulations.
#' - `tsim-submit.R`: 
#'     This script is a wrapper to submit jobs to the SLURM cluster.
#' 
#' The results of the simulations are saved as RDS files containing a list
#' of `SummarizedBenchmark` objects of length M.
#'
#' @author Patrick Kimes

## ##############################################################################
## parse simulation run parameters
## ##############################################################################

args <- commandArgs(trailingOnly = TRUE)

M <- as.integer(args[1])
ncores <- as.integer(args[2])
setting_vparam <- args[3]
setting_icparam <- args[4]


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
## devtools::install_github("areyesq89/SummarizedBenchmark")


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
source("R/du_tsim.R")

## parameter settings
source("tsim-core-settings.R")

## helper digest function
library("digest")

## set location
projdir <- "/n/irizarryfs01_backed_up/pkimes/project/project-02-benchmark-fdr"

## set parallelization
registerDoParallel(cores = ncores)


## ##############################################################################
## define baseline (null) simulation settings
## ##############################################################################

## setting 0: (base) null simulation setting
setting_base <- list(m = 20000,          # number of hypothesis tests
                     pi0 = 1,            # proportion of null hypotheses
                     effect_size = 0,    # expected mean diff of non-null tests
                     n_samples = 10,     # number of samples per group
                     n_groups = 2,       # number of groups in contrast
                     icovariate = FALSE) # self-explainatory


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


## ##############################################################################
## generate simulation settings
## ##############################################################################

settings <- define_settings(setting_base, setting_vparam, setting_icparam)


## ##############################################################################
## run and save simulations
## ##############################################################################

## save data and results if only 1 replication requested
if (M == 1) {

    for (idx in seq(settings)) {
        iset <- settings[[idx]]

        ## check if sim already run
        outf <- paste0("data/M", M, "/", setting_vparam, "-", setting_icparam, "/",
                       "results-", names(settings)[idx], "-M", M, ".rdata")
        if (file.exists(outf)) {
            next
        }
        dir.create(dirname(outf), showWarnings = FALSE, recursive = TRUE)
        
        sim_seed <- (as.integer(Sys.time()) %% 1e6)
        set.seed(sim_seed)

        ## simulate data with seed
        sim_df <- do.call(du_tsim, iset)
        names(sim_df)[which(names(sim_df) == "H")] <- "qvalue"
        
        ## calc data digest
        sim_digest <- digest::sha1(sim_df)
        
        ## create SummarizedBenchmark
        sb <- buildBench(bd, sim_df, truthCol = "qvalue", ptabular = TRUE)
        
        ## add simulation information to metadata
        metadata(sb)$sim_func <- du_tsim
        metadata(sb)$sim_parameters <- iset
        metadata(sb)$sim_seed <- sim_seed
        metadata(sb)$sim_digest <- sim_digest
        
        save(sb, sim_df, file = outf)
    }

} else { 

    ## loop over parameter settings
    for (idx in seq(settings)) {
        iset <- settings[[idx]]
        
        ## check if sim already run
        outf <- paste0("data/M", M, "/", setting_vparam, "-", setting_icparam, "/",
                       "results-", names(settings)[idx], "-M", M, ".rdata")
        if (file.exists(outf)) {
            next
        }
        dir.create(dirname(outf), showWarnings = FALSE, recursive = TRUE)
        
        sblist <- foreach(i = 1:M, .verbose = T) %dopar% {
            ## will break if i > ~2000 (integer seed value too large)
            sim_seed <- (as.integer(Sys.time()) %% 1e6) + (i * 1e6)
            set.seed(sim_seed)
            
            ## simulate data with seed
            sim_df <- do.call(du_tsim, iset)
            names(sim_df)[which(names(sim_df) == "H")] <- "qvalue"

            ## calc data digest
            sim_digest <- digest::sha1(sim_df)

            ## create SummarizedBenchmark
            sb <- buildBench(bd, sim_df, truthCol = "qvalue", ptabular = TRUE)

            ## add simulation information to metadata
            metadata(sb)$sim_func <- du_tsim
            metadata(sb)$sim_parameters <- iset
            metadata(sb)$sim_seed <- sim_seed
            metadata(sb)$sim_digest <- sim_digest

            sb
        }
        saveRDS(sblist, file = outf)
    }
}
