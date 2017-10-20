#' Simulation Setting Blocks for FDR Benchmarking
#'
#' Function takes a baseline set of simulation settings for the
#' `du_ttest_sim` function and returns a list of alternative simulation
#' settings which should be run jointly. 
#'
#' @param sbase
#' @param vparam variable parameter; should be one of: "pi0",
#'        "esize_fixed", "esize_random"
#' @param icparam independent/informative covariate parameter; should be
#'        one of: "uniform", "bl"
#'
#' @description
#' The following alternative parameter groups are defined.
#' - `sidx = 1`: no informative covariate, varying null proportion
#' 
#' @author Patrick Kimes

define_settings <- function(sbase, vparam, icparam) {

    ## check that variable parameter is one of specified set
    vp <- c("pi0", "esize_fixed", "esize_random_ua", "esize_random_shift", "altnoise")
    stopifnot(vparam %in% vp)
    
    ## check that indep/inform covariate parameter is one of specified set
    ic <- c("uniform", "bl")
    stopifnot(icparam %in% ic)

    ## filter out unsupported setting pairs
    if (vparam == "pi0" && icparam == "bl") {
        stop("cant run 'pi0' simulations w/ 'bl' format indep covariate")
    }

    ## ##########################################################################
    ## define type of informative/independent covariate by 'icparam'

    if (icparam == "uniform") {
        ## don't need to do anything
    } else if (icparam == "bl") {
        source("R/funcs_pi0.R")
        sbase <- replace(sbase, c("pi0", "icovariate"), c(pi0_smooth1, runif))
    }
    
    ## ##########################################################################
    ## define main simulations settings by 'vparam'

    if (vparam == "pi0") {
        ## varying pi0 (null proportions) #######################################
        settings <- lapply(seq(.1, 1, by=.1),
                           function(p) { replace(sbase,
                                                 c("pi0", "tstat"),
                                                 c(p, function(zz) { 1.5 }))
                           })
        names(settings) <- paste0("altp0_", 1:10)

    } else if (vparam == "esize_fixed") {
        ## varying effect size (fixed value) ####################################
        settings <- lapply(seq(0, 2, by=.2),
                           function(x) { replace(sbase,
                                                 c("pi0", "tstat"),
                                                 c(0.8, function(zz) { x })
                           })
        names(settings) <- paste0("alteff_", seq(0, 20, by=2))
        
    } else if (vparam == "esize_random_ua") {
        ## varying effect size (stochastic/UA) ##################################
        source("R/funcs_tstat.R")
        settings <- list("alteff_spiky" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_spiky)),
                         "alteff_skew" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_skew)),
                         "alteff_bignormal" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_bignormal)),
                         "alteff_bimodal" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_bimodal)),
                         "alteff_flat_top" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_flat_top)),
                         "alteff_near_normal" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, sampler_near_normal)))
        
    } else if (vparam == "esize_random_shift") {
        ## varying effect size not  ##################################
        ## add settings for rnorm/rt/rchisq_generator
        settings <- list("alteff_normal_shift1" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rnorm_generator(1))),
                         "alteff_normal_shift2" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rnorm_generator(2))),
                         "alteff_t_shift1_df10" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rt_generator(10, 1))),
                         "alteff_t_shift2_df10" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rt_generator(10, 1))),
                         "alteff_chisq_shift1_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 1))),
                         "alteff_chisq_shift2_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 1))))
        
    } else if (vparam == "altnoise") {
        ## add settings for rnorm/rt/rchisq_generator
        settings <- list("altnoise_t_df5" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, sampler_bimodal, rt_perturber(5), rt_2pvaluer(5))),
                         "altnoise_t_df10" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, sampler_bimodal, rt_perturber(10), rt_2pvaluer(10))),
                         "altnoise_chisq_df1" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 1+3^2), rchisq_perturber(1),
                                       rchisq_pvaluer(1))),
                         "altnoise_chisq_df4" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 1+3^2), rchisq_perturber(3),
                                       rchisq_pvaluer(1))))
    }

    return(settings)
}


