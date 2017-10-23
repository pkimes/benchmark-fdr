## ##############################################################################
## test all settings
## ##############################################################################
source("../R/du_psim.R")
source("../R/funcs_pi0.R")
source("../R/funcs_tstat.R")


## 1. pi0 fixed
dat <- du_psim(m = 1e4, pi0 = 0.8,
               tstat = rnorm_generator(2, 1),
               tstat_dist = rnorm_perturber(1),
               null_dist = rnorm_2pvaluer(1),
               icovariate = runif)

## 2. pi0 functional in informative covariate
dat <- du_psim(m = 1e4, pi0 = pi0_smooth1,
               tstat = rnorm_generator(2, 1),
               tstat_dist = rnorm_perturber(1),
               null_dist = rnorm_2pvaluer(1),
               icovariate = runif)

## 3. using t distributions
dat <- du_psim(m = 1e4, pi0 = pi0_smooth1,
               tstat = rnorm_generator(2, 1),
               tstat_dist = rt_perturber(3),
               null_dist = rt_2pvaluer(3),
               icovariate = runif)

## 4. using chisq distributions
dat <- du_psim(m = 1e4, pi0 = pi0_smooth1,
               tstat = rchisq_generator(1, 3^2),
               tstat_dist = rchisq_perturber(1), 
               null_dist = rchisq_pvaluer(3),
               icovariate = runif)




