#' Direct Parameter Simulation Settings for FDR Benchmarking
#'
#' Function takes a baseline set of simulation settings for the
#' `du_psim` function and returns a list of alternative simulation
#' settings which should be run jointly. 
#'
#' @param sbase default simulation settings to be modified.
#' @param vparam variable parameter; should be one of: "pi0",
#'        "esize_fixed", "esize_random_ua", "esize_random_shift",
#'        "altnoise", "allnull".
#' @param icparam independent/informative covariate parameter; should be
#'        one of: "uniform", "bl".
#'
#' @description
#' The following simulation groups are defined:
#' - `vparam`: variable parameter
#'     - `pi0`: for fixed effect size, varying null proportion.
#'     - `esize_fixed`: varying effect size taken to be the same fixed
#'                      value across all alternative tests.
#'     - `esize_random_ua`: varying effect size sampled from a (mostly)
#'                          unimodal set of distributions centered at zero
#'                          (distributions from ASH paper).
#'     - `esize_random_shift`: varying effect size sampled from a set of
#'                             distributions not centered at zero
#'                             (shifted normal, non-central t, non-central chi-sq). 
#'     - `altnoise`: varying sampling distributions for perturbing the
#'                   test statistics (both null and alternative)
#'                   (shifted normal, non-central t, non-central chi-sq).
#'     - `allnull`: similar to `altnoise` setting, but with pi0 = 1, i.e. all test
#'                  statistics are simulated under the null distribution.
#' - `icparam`: independent covariate structure
#'     - `uniform`: covariate sampled from uniform (0, 1) interval AND
#'                  no association between covariate and distribution of
#'                  p-values is introduced.
#'     - `bl`: covariate sampled from uniform (0, 1) interval AND
#'             the probability of a test being null (pi0) is given a functional
#'             form of the single covariate (as in the Boca-Leek paper).
#' 
#' @author Patrick Kimes

psim_settings <- function(sbase, vparam, icparam) {

    ## check that variable parameter is one of specified set
    vp <- c("pi0", "esize_fixed", "esize_random_ua", "esize_random_shift",
            "altnoise", "allnull")
    stopifnot(vparam %in% vp)
    
    ## check that indep/inform covariate parameter is one of specified set
    ic <- c("uniform", "bl")
    stopifnot(icparam %in% ic)

    ## filter out unsupported setting pairs
    if (vparam == "pi0" && icparam == "bl") {
        stop("cant run 'pi0' simulations w/ 'bl' format indep covariate")
    }
    if (vparam == "allnull" && icparam == "bl") {
        stop("cant run 'allnull' simulations w/ 'bl' format indep covariate")
    }
    
    ## ##########################################################################
    ## define main simulations settings by 'vparam'
    
    if (vparam == "pi0") {
        ## varying pi0 (null proportions) #######################################
        settings <- lapply(seq(.1, 1, by=.1),
                           function(p) { replace(sbase,
                                                 c("pi0", "tstat"),
                                                 c(p, function(zz) { 3.0 }))
                           })
        names(settings) <- paste0("altp0_", 1:10)
        
    } else if (vparam == "esize_fixed") {
        ## varying effect size (fixed value) ####################################
        settings <- lapply(seq(0, 5, by=.5),
                           function(x) { replace(sbase,
                                                 c("pi0", "tstat"),
                                                 c(0.8, function(zz) { x }))
                           })
        names(settings) <- paste0("alteff_", seq(0, 50, by=5))
        
    } else if (vparam == "esize_random_ua") {
        ## varying effect size (stochastic/UA) ##################################
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
        ## varying effect size (mean shifted distributions) #####################
        settings <- list("alteff_normal_shift1" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rnorm_generator(1))),
                         "alteff_normal_shift2" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rnorm_generator(2))),
                         "alteff_normal_shift3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rnorm_generator(3))),
                         "alteff_t_shift1_df10" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rt_generator(10, 1))),
                         "alteff_t_shift2_df10" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rt_generator(10, 2))),
                         "alteff_t_shift3_df10" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rt_generator(10, 3))),
                         "alteff_chisq_shift0_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 0))),
                         "alteff_chisq_shift1_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 1))),
                         "alteff_chisq_shift2_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 2))),
                         "alteff_chisq_shift3_df3" =
                             replace(sbase, c("pi0", "tstat"),
                                     c(0.8, rchisq_generator(3, 3))))
        
    } else if (vparam == "altnoise") {
        ## varying sampling noise ###############################################
        settings <- list("altnoise_t_df5_bimodal" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, sampler_bimodal, rt_perturber(5), rt_2pvaluer(5))),
                         "altnoise_t_df10_bimodal" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, sampler_bimodal, rt_perturber(10), rt_2pvaluer(10))),
                         "altnoise_chisq_df1_shift2sq" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 2^2), rchisq_perturber(1),
                                       rchisq_pvaluer(1))),
                         "altnoise_chisq_df4_shift2sq" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 2^2), rchisq_perturber(4),
                                       rchisq_pvaluer(4))),
                         "altnoise_chisq_df1_shift3sq" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 3^2), rchisq_perturber(1),
                                       rchisq_pvaluer(1))),
                         "altnoise_chisq_df4_shift3sq" =
                             replace(sbase, c("pi0", "tstat", "tstat_dist", "null_dist"),
                                     c(0.8, rchisq_generator(1, 3^2), rchisq_perturber(4),
                                       rchisq_pvaluer(4))))

            
    } else if (vparam == "allnull") {
        ## varying sampling noise ###############################################
        settings <- list("altnoise_normal_null" =
                             replace(sbase, c("tstat", "tstat_dist", "null_dist"),
                                     c(function(x) { 0 }, rnorm_perturber(1), rnorm_2pvaluer(1))),
                         "altnoise_t_df5_null" =
                             replace(sbase, c("tstat", "tstat_dist", "null_dist"),
                                     c(function(x) { 0 }, rt_perturber(5), rt_2pvaluer(5))),
                         "altnoise_t_df10_null" =
                             replace(sbase, c("tstat", "tstat_dist", "null_dist"),
                                     c(function(x) { 0 }, rt_perturber(10), rt_2pvaluer(10))),
                         "altnoise_chisq_df1_null" =
                             replace(sbase, c("tstat", "tstat_dist", "null_dist"),
                                     c(function(x) { 0 }, rchisq_perturber(1),
                                       rchisq_pvaluer(1))),
                         "altnoise_chisq_df4_null" =
                             replace(sbase, c("tstat", "tstat_dist", "null_dist"),
                                     c(function(x) { 0 }, rchisq_perturber(4),
                                       rchisq_pvaluer(4))))
    }
    
    ## ##########################################################################
    ## define type of informative/independent covariate by 'icparam'
    
    if (icparam == "uniform") {
        ## don't need to do anything
    } else if (icparam == "bl") {
        settings <- lapply(settings, replace, c("pi0", "icovariate"), c(pi0_smooth1, runif))
    }
    
    return(settings)
}


