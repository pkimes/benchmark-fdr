#' Simulate Collection of T-Tests
#'
#' Function for simulating a large set of (independent) t-tests. Simulation is performed
#' at the "raw data" level. Test results are returned as a data.frame with one row per
#' test. 
#' 
#' @param m number of hypothesis tests. 
#' @param pi0 proportion of true null hypotheses, can either be a single numeric value
#'        between 0 and 1 or a function which takes a numeric vector and returns a
#'        vector of the same length containing null probabilities between 0 and 1. 
#' @param effect_size expected mean difference, can either be a single numeric value
#'        or a function which takes an integer as input and returns a numeric vector
#'        of the specified length.
#' @param n_samples integer vector specifying the number of samples per group.
#'        If `length(n_samples) == 1`, then the same number of samples is used
#'        for the `n_groups` groups. (default = 10)
#' @param n_groups number of groups in comparison. Currently, only `n_groups = 2`
#'        is supported. (default = 2)
#' @param icovariate specification for the independent covariate, can either be a logical
#'        value or a function which takes an integer value and returns a vector of the
#'        specified length to be used as the independent covariate. If logical and TRUE, 
#'        then the variance of the pooled samples is returned as the independent covariate 
#'        as used in the IHW simulations, and if logical and FALSE, the an uniformative
#'        covariate (randomly sample from Uniform[0, 1]) is returned as the independent
#'        covariate. (default = TRUE)
#' @param seed integer seed for random number generator, ignored if NULL (default = NULL) 
#'
#' @return
#' data.frame of test results for `m` simulated data sets, with columns:
#' * `H`: 0/1 indicator whether data simulated under null (0) or alternative (1)
#' * `test_statistic`: t-test statistic
#' * `effect_size`: group mean difference
#' * `pval`: t-test p-value
#' * `ind_covariate`: the independent covariate 
#' * `SE`: standard error across all samples (estimated using the ratio of difference in means / t-test statistic)
#'
#' @details
#' If a function is specified for either of `effect_size` or `icovariate`, a function must also be
#' specified for the other parameter.
#' This function builds on `IHWpaper::du_ttest_sim` originally written by Nikos Ignatiadis
#' and released under the Artistic-2.0 license. The original function is available at:
#' https://github.com/nignatiadis/IHWpaper
#'
#' @md
#' @importFrom genefilter rowtests
#' @importFrom genefilter rowSds
#' @author Stephanie Hicks, Patrick Kimes
du_tsim <- function(m, pi0, effect_size, n_samples = 10, n_groups = 2, 
                    icovariate = TRUE, seed = NULL) {
    ## use specified random seed
    if (!is.null(seed)) {
        set.seed(seed)
    }

    ## check specified number of groups and group sizes
    if (length(n_samples) == 1) {
        n_samples <- rep(n_samples, n_groups)
    }
    stopifnot(n_groups == 2)
    stopifnot(length(n_samples) == n_groups)
    ntotal <- sum(n_samples)

    ## check for functional icovariate, pi0 specification
    if (is.function(icovariate) && is.function(pi0)) {
        ## simulate indep covariate from icovariate function
        ind_cov <- icovariate(m)
        ## pi0 returns probability of null, sample alts from [1 - pi0s]
        pi0s <- pi0(ind_cov)
        stopifnot(length(pi0s) == m)
        alts <- rbinom(m, 1, 1 - pi0s)
    } else if (is.logical(icovariate) && is.numeric(pi0)) {
        ## use flat pi0 rate
        m0 <- ceiling(m * pi0)
        alts <- sample(1:m, m - m0)
    } else {
        stop("icovariate and pi0 must either both be functions or ",
             "logical and numeric, respectively.")
    }

    ## check for functional effect_size specification
    if (is.function(effect_size)) {
        esize <- effect_size
    } else if (is.numeric(effect_size) && length(effect_size) == 1) {
        esize <- rep(effect_size, length(alts))
    } else {
        stop("effect_size must either be a function or ",
             "a numeric vector of length 1")
    }
    
    ## simulate data from mixture normal
    z_table <- matrix(rnorm(n_total * m), ncol = n_total)
    z_table[alts, (n_samples / 2 + 1):n_samples] <-
        z_table[alts, (n_samples / 2 + 1):n_samples] + effect_size
    
    ## if ind_cov not simulated, use sample cov or uniform (as in IHW)
    if (is.logical(icovariate)) {
        if(icovariate) {
            ind_cov <- genefilter::rowSds(z_table)
        } else {
            ind_cov <- runif(m, 0, 1)
        }
    }
    
    ## null/alt indicator
    H <- rep(0, m)
    H[alts] <- 1

    ## perform t-tests
    gF <- factor(rep(1:n_groups, each = n_samples / n_groups))
    t_test <- genefilter::rowttests(z_table, gF)

    ## standard error needed for ASH
    SE <- t_test$dm / t_test$statistic

    ## return test results as data.frame
    data.frame(H = H, test_statistic = t_test$statistic, effect_size = t_test$dm, 
               pval = t_test$p.value, ind_covariate = sds, SE = SE)
}
