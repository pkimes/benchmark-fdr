## ##############################################################################
## test all settings
## ##############################################################################
source("du_psim.R")
source("funcs_pi0.R")
source("funcs_tstat.R")


## 1. pi0 fixed
dat <- du_psim(m = 1e4, pi0 = 0.8,
               tstat = rnorm_generator(2, 1),
               tstat_dist = rnorm_perturber(1),
               null_dist = rnorm_2pvaluer(1),
               icovariate = runif)

## 1. pi0 functional in informative covariate
dat <- du_psim(m = 1e4, pi0 = pi0_smooth1,
               tstat = rnorm_generator(2, 1),
               tstat_dist = rnorm_perturber(1),
               null_dist = rnorm_2pvaluer(1),
               icovariate = runif)


