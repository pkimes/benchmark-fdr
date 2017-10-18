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
rnormmix_sampler <- function(g) {
    k <- unique(sapply(g, length))
    stopifnot(length(k) == 1)
    stopifnot(sum(g[[1]]) == 1)
    function(n) {
        comp <- sample(1:k, n, g[[1]], replace=TRUE)
        rnorm(n, g[[2]][comp], g[[3]][comp])
    }
}

## ASH effect-size distributions
sampler_spiky <- rnormmix_sampler(params$spiky)
sampler_skew <- rnormmix_sampler(params$skew)
sampler_bignormal <- rnormmix_sampler(params$bignormal)
sampler_bimodal <- rnormmix_sampler(params$bimodal)
sampler_flat_top <- rnormmix_sampler(params$flat_top)
sampler_near_normal <- rnormmix_sampler(params$near_normal)

## ##############################################################################
## END CODE DERIVED FROM EXTERNAL SOURCE
## ##############################################################################
