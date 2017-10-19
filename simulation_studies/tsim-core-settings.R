#' Simulation Setting Blocks for FDR Benchmarking
#'
#' Function takes a baseline set of simulation settings for the
#' `du_ttest_sim` function and returns a list of alternative simulation
#' settings which should be run jointly. 
#'
#' @param sbase
#' @param vparam variable parameter; should be one of: "pi0", "n",
#'        "esize_fixed", "esize_random"
#' @param icparam independent/informative covariate parameter; should be
#'        one of: "uniform", "se", "bl"
#'
#' @description
#' The following alternative parameter groups are defined.
#' - `sidx = 1`: no informative covariate, varying null proportion
#' 
#' @author Patrick Kimes

define_settings <- function(sbase, vparam, icparam) {

    ## check that variable parameter is one of specified set
    vp <- c("pi0", "n", "esize_fixed", "esize_random")
    stopifnot(vparam %in% vp)
    
    ## check that indep/inform covariate parameter is one of specified set
    ic <- c("uniform", "se", "bl")
    stopifnot(icparam %in% ic)

    ## filter out unsupported setting pairs
    if (vparam == "pi0" && icparam == "bl") {
        stop("cant run 'pi0' simulations w/ 'bl' format indep covariate")
    }

    ## ##########################################################################
    ## define type of informative/independent covariate by 'icparam'

    if (icparam == "uniform") {
        ## don't need to do anything
    } else if (icparam == "se") {
        sbase <- replace(sbase, "icovariate", TRUE)
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
                                                 c("pi0", "effect_size"),
                                                 c(p, 1.5))
                           })
        names(settings) <- paste0("altp0_", 1:10)

    } else if (vparam == "esize_fixed") {
        ## varying effect size (fixed value) ####################################
        settings <- lapply(seq(0, 2, by=.2),
                           function(x) { replace(sbase,
                                                 c("pi0", "effect_size"),
                                                 c(0.8, x))
                           })
        names(settings) <- paste0("alteff_", seq(0, 20, by=2))

    } else if (vparam == "n") {
        ## varying sample size ##################################################
        settings <- list("altn_3" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(3, 0.8, 1.5)),
                         "altn_6" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(6, 0.8, 1.5)),
                         "altn_20" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(20, 0.8, 1.5)),
                         "altn_6_12" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(c(6, 12), 0.8, 1.5))
                         )
        
    } else if (vparam == "esize_random") {
        ## varying effect size (stochastic/UA) ##################################
        source("R/funcs_tstat.R")
        settings <- list("alteff_spiky" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_spiky)),
                         "alteff_skew" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_skew)),
                         "alteff_bignormal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_bignormal)),
                         "alteff_bimodal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_bimodal)),
                         "alteff_flat_top" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_flat_top)),
                         "alteff_near_normal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(0.8, sampler_near_normal))
                         )
    }
        
    return(settings)
}


