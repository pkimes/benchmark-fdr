#' Simulate Collection of Test Statistics
#'
#' Function for simulating a large set of (independent) test results.
#' Unlike `du_tsim.R` which simulates actual data sets to which two-sample
#' t-tests are applied, this function directly simulates test statistics
#' and p-values without any underlying data.
#' This simulator more closely represents the set of simulations carried
#' out in the manuscripts of the Boca-Leek, ASH, and IHW methods.
#' (For IHW, in particular, the simulations illustrating the "size investing strategy"
#' fall under this class of simulations.)
#' 
#' @param m number of hypothesis tests. 
#' @param pi0 proportion of true null hypotheses, can either be a single numeric value
#'        between 0 and 1 or a function which takes a numeric vector and returns a
#'        vector of the same length containing null probabilities between 0 and 1. 
#' @param icovariate specification for the independent covariate, must be a function
#'        which takes an integer value and returns a vector of the specified length to
#'        be used as the independent covariate.
#' @param tstat expected test statistic under the alternative, can either be a single
#'        numeric value or a function which takes an integer as input and returns a
#'        numeric vector of test statistics of the specified length.
#' @param tstat_dist sampling distribution of test statistic, must be a function
#'        which takes a vector of test statistics, and returns a list of perturbed
#'        test statistics of the same length.
#' @param null_dist distribution under the null used to calculate p-values from the
#'        test statistic, must be a function which takes the full vector of null and
#'        alternative test statistics and return the corresponding p-value. 
#' @param seed integer seed for random number generator, ignored if NULL. (default = NULL) 
#'
#' @return
#' data.frame of test results for `m` simulated data sets, with columns:
#' * `H`: 0/1 indicator whether data simulated under null (0) or alternative (1)
#' * `test_statistic`: t-test statistic
#' * `effect_size`: group mean difference
#' * `pval`: t-test p-value
#' * `ind_covariate`: the independent covariate 
#'
#' @details
#' If a function is specified for either of `effect_size` or `icovariate`, a function must also be
#' specified for the other parameter.
#'
#' @md
#' @author Patrick Kimes
du_psim <- function(m, pi0, tstat, tstat_dist, null_dist, icovariate, seed = NULL) {

    library("genefilter")
    
    ## use specified random seed
    if (!is.null(seed)) {
        set.seed(seed)
    }

    stopifnot(is.function(tstat))
    stopifnot(is.function(tstat_dist))
    stopifnot(is.function(icovariate))
              
    ## convert pi0 to function
    if (!is.function(pi0)) {
        pi0 <- function(x) { pi0 }
    }
    
    ## simulate indep covariate from icovariate function
    ind_cov <- icovariate(m)
    stopifnot(length(ind_cov) == m)

    ## pi0 returns probability of null, sample alts from [1 - pi0s]
    if (is.function(pi0)) {
        pi0s <- pi0(ind_cov)
    } else if (length(pi0) == 1) {
        pi0s <- rep(pi0, m)
    } else {
        stop("pi0 must be function or single numeric value")
    }
    stopifnot(length(pi0s) == m)
    stopifnot(min(pi0s) >= 0 && max(pi0s) <= 1)
    alts <- which(rbinom(m, 1, 1 - pi0s) == 1)
    
    ## generate set of null (0) and alternative test statistics
    ts <- rep(0, m)
    ts[alts] <- tstat(length(alts))
    ts <- tstat_dist(ts)
    stopifnot(length(ts) == length(alts))
    
    ## null/alt indicator
    H <- rep(0, m)
    H[alts] <- 1

    ## calculate p-values and 'effect size' (just N(0,1) z-score)
    pv <- null_dist(ts)
    es <- qnorm(pv)
    
    ## return test results as data.frame
    data.frame(H = H, test_statistic = ts, effect_size = es,
               pval = pv, ind_covariate = ind_cov)
}
