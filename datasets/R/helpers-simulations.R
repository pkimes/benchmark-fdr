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
#' * `test_statistic`: simulated test statistic (NOT scaled by SE estimate)
#' * `effect_size`: z-score transform of p-value (not really effect size)
#' * `pval`: test p-value calculated from test-statistic using `null_dist`
#' * `ind_covariate`: the independent covariate
#' * `SE`: true standard deviations for sampling distributions (for ASH)
#'
#' @details
#' If a function is specified for either of `effect_size` or `icovariate`, a function must also be
#' specified for the other parameter.
#'
#' @md
#' @author Patrick Kimes
du_psim <- function(m, pi0, tstat, tstat_dist, null_dist, icovariate, seed = NULL) {
    
    ## use specified random seed
    if (!is.null(seed)) {
        set.seed(seed)
    }

    stopifnot(is.function(tstat))
    stopifnot(is.function(tstat_dist))
    stopifnot(is.function(icovariate))
    
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
    if (length(alts) > 0) {
        ts[alts] <- tstat(length(alts))
    }
    
    ## determine (true) SD of sampling dist used to perturb each stat
    SE <- tstat_dist(ts, se = TRUE)

    ## perturb each stat
    ts <- tstat_dist(ts)
    stopifnot(length(ts) == m)
    
    ## null/alt indicator
    H <- rep(0, m)
    H[alts] <- 1

    ## calculate p-values and 'effect size' (just N(0,1) z-score)
    pv <- null_dist(ts)
    es <- qnorm(1 - pv / 2) * sign(ts)

    ## return test results as data.frame
    data.frame(H = H, test_statistic = ts, effect_size = es,
               pval = pv, ind_covariate = ind_cov, SE = SE)
}


## ##############################################################################
## START CODE FROM EXTERNAL SOURCE
## source: https://github.com/SiminaB/Fdr-regression/blob/71b123f/functions.R
## license: not listed
## date: 2017/10/17
## modified: Patrick Kimes
## ##############################################################################

##------Functions of covariates-------##
pi0_f1 <- function(x) {
    p2 <- -0.2
    p1 <- 1.2
    a <- 4/(p1-p2)^2
    
    y <- -a*(x-p1)*(x-p2)
    y[x >= 0.7] <- -a*(0.7-p1)*(0.7-p2)
    y[x <= (p1+p2)/2] <- 1
    y
}
pi0_f2 <- function(x) {
    y <- rep(0, length=length(x))
    y[x >= 0.7] <- -2.5*(x[x >= 0.7]-0.7)^2
    y
}
pi0_f3 <- function(x) {
    y <- rep(0, length=length(x))
    y[x < 0.7] <- -(x[x < 0.7]-0.1)^2
    y[x >= 0.7] <- -(min(x[x >= 0.7])-0.1)^2
    y[x<=0.1] <- 0
    y
}

##smooth function of one covariate for different levels of second covariate
pi0_smooth2 <- function(x1, x2) {
    y1 <- pi0_f1(x1)
    y2 <- pi0_f2(x1)
    y3 <- pi0_f3(x1)
    y <- rep(0, length(x1))
    y[x2 == 1] <- y1[x2 == 1] + y2[x2 == 1] + 0.12*y3[x2 == 1]
    y[x2 == 2] <- y1[x2 == 2] + 0.5*y2[x2 == 2] + 0.06*y3[x2 == 2]
    y[x2 == 3] <- y1[x2 == 3] + 0.3*y2[x2 == 3] 
    y
}

##smooth function of a single variable
##note: assumes x in (0, 1)
pi0_smooth1 <- function(x) {
    y1 <- pi0_f1(x)
    y2 <- pi0_f2(x)
    y3 <- pi0_f3(x)
    y <- rep(0, length(x))
    y <- y1 + y2 + 0.12*y3
    y
}

## ##############################################################################
## END CODE FROM EXTERNAL SOURCE
## ##############################################################################

## ##############################################################################
## alternative Boca-Leek style informative covariates
## - if x ~ U(0,1), then expect 80%, 90%, 95% pi0 for all cases
## ##############################################################################

## step function (4 steps)
pi0_step <- function(x, pi0) {
    stopifnot(pi0 < 1, pi0 >= 0.5)
    pi0 - (1-pi0)/2 +
        (1-pi0)/4 * (x > 0.25) +
        (1-pi0)/2 * (x > 0.5) +
        (1-pi0)/4 * (x > 0.75)
}

## shifted/stretched cubic function
pi0_cubic <- function(x, pi0) {
    stopifnot(pi0 < 1, pi0 >= 2/3)
    (1-x)^(1/3) * 4*(1-pi0) + 4 * pi0 - 3
}

## shifted/stretched cosine function (non-monotone)
pi0_cosine <- function(x, pi0) {
    stopifnot(pi0 < 1, pi0 >= 1/2)
    (1-pi0) * cos(2*pi*x) + pi0
}

pi0_sine <- function(x, pi0) {
    stopifnot(pi0 < 1, pi0 >= 1/2)
    (1-pi0) * sin(2*pi*x) + pi0
}


## ##############################################################################
## CODE DERIVED FROM EXTERNAL SOURCE
## source: https://github.com/stephenslab/ash/blob/94a3347/code/dsc-shrink/datamakers/datamaker.R
## source: https://github.com/stephenslab/ash/blob/d06047e/code/dsc-shrink/add_named_scenarios.R
## license: not listed
## date: 2017/10/18
## modified: Patrick Kimes
## ##############################################################################

## test statistic distributions used in ASH simulations
params <- list(spiky = list(c(.4, .2, .2, .2), c(0, 0, 0, 0), c(.25, .5, 1, 2)),
               skew = list(c(1/4, 1/4, 1/3, 1/6), c(-2, -1, 0, 1), c(2, 1.5, 1, 1)),
               bignormal = list(c(1), c(0), c(4)),
               bimodal = list(c(0.5, 0.5), c(-2, 2), c(1, 1)),
               flat_top = list(rep(1/7, 7), c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), rep(0.5, 7)),
               near_normal = list(c(2/3, 1/3), c(0, 0), c(1, 2)))

## function factory to generate mixture normal samplers 
rnormmix_generator <- function(g) {
    k <- unique(sapply(g, length))
    stopifnot(length(k) == 1)
    stopifnot(sum(g[[1]]) == 1)
    function(n) {
        comp <- sample(1:k, n, g[[1]], replace=TRUE)
        rnorm(n, g[[2]][comp], g[[3]][comp])
    }
}

## ASH effect-size distributions
sampler_spiky <- rnormmix_generator(params$spiky)
sampler_skew <- rnormmix_generator(params$skew)
sampler_bignormal <- rnormmix_generator(params$bignormal)
sampler_bimodal <- rnormmix_generator(params$bimodal)
sampler_flat_top <- rnormmix_generator(params$flat_top)
sampler_near_normal <- rnormmix_generator(params$near_normal)

## ##############################################################################
## END CODE DERIVED FROM EXTERNAL SOURCE
## ##############################################################################

## function factory: normal distribution
rnorm_generator <- function(m, s = 1) {
    function(n) {
        rnorm(n, m, s)
    }
}

## function factory: non-central t distribution
rt_generator <- function(df, m = 0) {
    if (m == 0) {
        return(function(n) { rt(n, df) })
    } else if (df > 1) {
        ncp <- m / gamma((df-1)/2) * gamma(df/2) / sqrt(df/2)
        return(function(n) { rt(n, df, ncp) })
    } else {
        stop("eep! my author didn't write me to accept ",
             "(df <= 1) & (m != 0). sorry!")
    }
}

## function factory: non-central chi-sq distribution
rchisq_generator <- function(df, ncp = 0) {
    return(function(n) { rchisq(n, df, ncp) })
}

## ##############################################################################
## ##############################################################################

## function factory: normal perturbation
## - returns functions which takes input vector and simulates normal
##   random variables with specified vector as means.
## - returned function can also return (true) SD of sampling dist
rnorm_perturber <- function(s = 1) {
    function(m, se = FALSE) {
        if (se) { return(rep(s, length(m))) }
        rnorm(length(m), m, s)
    }
}

## function factory: non-central t distribution
## - returns functions which takes input vector and simulates non-central t
##   random variables with specified vector as means.
## - returned function can also return (true) SD of sampling dist
rt_perturber <- function(df) {
    if (df <= 2) {
        stop("var not defined for non-central t with df <= 2.")
    }
    function(m, se = FALSE) {
        ncp <- m / gamma((df-1)/2) * gamma(df/2) / sqrt(df/2)
        if (se) { return(sqrt(df*(1+ncp^2) / (df-2) - m^2)) }
        ## rt(..) w/ length(ncp) > 1 throws warning - just call underlying code
        rnorm(length(m), ncp) / sqrt(rchisq(length(m), df) / df)
    }
}

## function factory: non-central chi-sq perturbation
## - returns functions which takes input vector and sim non-central chi-sq
##   random variables with specified vector as means.
## - returned function can also return (true) SD of sampling dist
rchisq_perturber <- function(df) {
    function(ncp, se = FALSE) {
        stopifnot(ncp >= 0)
        if (se) { return(sqrt(2 * df + 4 * ncp)) }
        rchisq(length(ncp), df, ncp)
    }
}

## ##############################################################################
## ##############################################################################

## function factory: two-sided normal p-value calculator
## @param s standard deviation of null zero-centered normal distribution
rnorm_2pvaluer <- function(s) {
    function(x) { 2 * (1 - pnorm(abs(x), 0, s)) }
}


## function factory: two-sided t p-value calculator
## @param df degrees of freedom of null zero-centered t-distribution
rt_2pvaluer <- function(df) {
    function(x) { 2 * (1 - pt(abs(x), df)) }
}


## function factory: chi-sq p-value calculator
## @param df degrees of freedom of null zero-centered t-distribution
rchisq_pvaluer <- function(df) {
    function(x) { 1 - pchisq(x, df) }
}

