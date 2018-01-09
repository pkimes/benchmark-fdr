#' T-Test Simulation Settings for FDR Benchmarking
#'
#' Function takes a baseline set of simulation settings for the
#' `du_tsim` function and returns a list of alternative simulation
#' settings which should be run jointly. 
#'
#' @param sbase default simulation settings to be modified.
#' @param vparam variable parameter; should be one of: "pi0", "n",
#'        "esize_fixed", "esize_random".
#' @param icparam independent/informative covariate parameter; should be
#'        one of: "uniform", "se", "bl".
#' @param base_pi0 pi0 value to use for simulations if `vparam` is one
#'        of: "esize_fixed", "esize_random_ua", "esize_random_shift",
#'        "altnoise". If `vparam` is either of `pi0` or `allnull`, then
#'        the value is ignored. (default = 0.8)
#'
#' @details
#' The following simulation groups are defined:
#' - `vparam`: variable parameter
#'     - `pi0`: for fixed effect size and sample size, varying null prop.
#'     - `n`: for fixed effect size and null prop, varying sample sizes.
#'     - `esize_fixed`: for fixed sample size and null prop, varying effect size
#'                      taken to be the same fixed value across all alternative tests.
#'     - `esize_random`: for fixed sample size and null prop, varying effect size
#'                       sampled from (mostly) unimodal set of distributions centered
#'                       centered at zero (distributions from ASH paper).
#' - `icparam`: independent covariate structure
#'     - `uniform`: covariate sampled from uniform (0, 1) interval.
#'     - `se`: standard error used as covariate (as in IHW paper).
#'     - `bl`: covariate sampled from uniform (0, 1) interval AND
#'             probability of test being null (pi0) takes a functional
#'             form of the single covariate (as in the Boca-Leek paper).
#'     - `bl-step-less`: covariate sampled from uniform (0, 1) interval AND
#'                       the probability of a test being null (pi0) is given a
#'                       functional form of the single covariate.
#'                       The functional form is a step function taking values in
#'                       {0.70, 0.75, 0.85, 0.90}. The step function is chosen
#'                       so the marginal pi0 is 0.80. 
#'     - `bl-step-more`: covariate sampled from uniform (0, 1) interval AND
#'                       the probability of a test being null (pi0) is given a
#'                       functional form of the single covariate.
#'                       The functional form is a step function taking values in
#'                       {0.6, 0.7, 0.9, 1.0}. The step function is chosen
#'                       so the marginal pi0 is 0.80.
#'     - `bl-cubic`: covariate sampled from uniform (0, 1) interval AND
#'                   the probability of a test being null (pi0) is given a
#'                   functional form of the single covariate.
#'                   The functional form is a reflected, stretched, and shifted
#'                   cubic function of the single covariate. The transformation
#'                   is chosen so the marginal pi0 is 0.80.
#'     - `bl-step-90`: same as `bl-step-less` but marginal pi0 of 0.90.
#'     - `bl-step-95`: same as `bl-step-less` but marginal pi0 of 0.95.
#'     - `bl-cubic-90`: same as `bl-cubic` but marginal pi0 of 0.90.
#'     - `bl-cubic-95`: same as `bl-cubic` but marginal pi0 of 0.95.
#' 
#' @author Patrick Kimes

tsim_settings <- function(sbase, vparam, icparam, base_pi0 = 0.8) {

    ## check that variable parameter is one of specified set
    vp <- c("pi0", "n", "esize_fixed", "esize_random")
    stopifnot(vparam %in% vp)
    
    ## check that indep/inform covariate parameter is one of specified set
    ic <- c("uniform", "se", "bl", "bl-step-less", "bl-step-more", "bl-cubic",
            "bl-step-90", "bl-step-95", "bl-cubic-90", "bl-cubic-95")
    stopifnot(icparam %in% ic)

    ## filter out unsupported setting pairs
    if (vparam == "pi0" && grepl("^bl", icparam)) {
        stop("cant run 'pi0' simulations w/ 'bl' format indep covariate")
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
                                                 c(base_pi0, x))
                           })
        names(settings) <- paste0("alteff_", seq(0, 20, by=2))

    } else if (vparam == "n") {
        ## varying sample size ##################################################
        settings <- list("altn_3" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(3, base_pi0, 1.5)),
                         "altn_6" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(6, base_pi0, 1.5)),
                         "altn_20" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(20, base_pi0, 1.5)),
                         "altn_6_12" =
                             replace(sbase, c("n_samples", "pi0", "effect_size"),
                                     c(c(6, 12), base_pi0, 1.5))
                         )
        
    } else if (vparam == "esize_random") {
        ## varying effect size (stochastic/UA) ##################################
        settings <- list("alteff_spiky" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_spiky)),
                         "alteff_skew" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_skew)),
                         "alteff_bignormal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_bignormal)),
                         "alteff_bimodal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_bimodal)),
                         "alteff_flat_top" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_flat_top)),
                         "alteff_near_normal" =
                             replace(sbase, c("pi0", "effect_size"),
                                     c(base_pi0, sampler_near_normal))
                         )
    }
        

    ## ##########################################################################
    ## define type of informative/independent covariate by 'icparam'

    if (icparam == "uniform") {
        ## don't need to do anything
    } else if (icparam == "se") {
        settings <- lapply(settings, replace, "icovariate", TRUE)
    } else if (icparam == "bl") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_smooth1, runif))
    } else if (icparam == "bl-step-less") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_step_less, runif))
    } else if (icparam == "bl-step-more") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_step_more, runif))
    } else if (icparam == "bl-cubic") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_cubic, runif))

    } else if (icparam == "bl-step-90") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_step90, runif))
    } else if (icparam == "bl-step-95") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_step95, runif))
    } else if (icparam == "bl-cubic-90") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_cubic90, runif))
    } else if (icparam == "bl-cubic-95") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_cubic95, runif))
    } 
    
    return(settings)
}


