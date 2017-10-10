#' Simulation Setting Blocks for FDR Benchmarking
#'
#' Function takes a baseline set of simulation settings for the
#' `du_ttest_sim` function and returns a list of alternative simulation
#' settings which should be run jointly. 
#'
#' @param sbase
#' @param sidx alternative setting index
#'
#' @description
#' The following alternative parameter groups are defined.
#' - `sidx = 1`: no informative covariate, varying null proportion
#' 
#' @author Patrick Kimes

define_settings <- function(sbase, sidx) {

    if (sidx == 1) { ## 1: varying pi0 (null proportions)
        settings <- lapply(seq(.1, 1, by=.1),
                           function(p) { replace(sbase,
                                                 c("pi0", "effect_size"),
                                                 c(p, 1.5))
                           })
        names(settings) <- paste0("altp0_", 1:10)

    } else if (sidx == 2) { ## 2: varying effect size
        settings <- lapply(seq(0, 2, by=.2),
                           function(x) { replace(sbase,
                                                 c("pi0", "effect_size"),
                                                 c(0.8, x))
                           })
        names(settings) <- paste0("alteff_", seq(0, 20, by=2))
        
    } else { ## not valid
        stop("Specified sidx is not valid!")
    }
    
    return(settings)
}
