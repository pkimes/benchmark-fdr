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
