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
rchisq_generator <- function(df, m = 0) {
    if (m == 0) {
        return(function(n) { rchisq(n, df) })
    } else {
        ncp <- m - df
        stopifnot(ncp > 0)
        return(function(n) { rchisq(n, df, ncp) })
    }
}

## ##############################################################################
## ##############################################################################

## function factory: normal perturbation
## - returns functions which takes input vector and simulates normal
##   random variables with specified vector as means.
rnorm_perturber <- function(s = 1) {
    function(m) {
        rnorm(length(m), m, s)
    }
}

## function factory: non-central t distribution
## - returns functions which takes input vector and simulates non-central t
##   random variables with specified vector as means.
rt_perturber <- function(df) {
    if (df <= 1) {
        stop("eep! mean for non-central t-distribution not defined for df = 1.")
    }
    
    function(m) {
        ncp <- m / gamma((df-1)/2) * gamma(df/2) / sqrt(df/2)
        rt(length(m), df, ncp)
    }
}

## function factory: non-central chi-sq perturbation
## - returns functions which takes input vector and sim non-central chi-sq
##   random variables with specified vector as means.
rchisq_perturber <- function(df) {
    function(m) {
        ncp <- m - df
        stopifnot(ncp > 0)
        rchisq(length(m), df, ncp)
    }
}

## ##############################################################################
## ##############################################################################

## function factory: two-sided normal p-value calculator
rnorm_2pvaluer <- function(s) {
    function(x) { 2 * (1 - pnorm(abs(x), 0, s)) }
}


## function factory: two-sided t p-value calculator
rt_2pvaluer <- function(df) {
    function(x) { 2 * (1 - pt(abs(x), df)) }
}


## function factory: chi-sq p-value calculator
rchisq_pvaluer <- function(df) {
    function(x) { 1 - pchisq(x, df) }
}

