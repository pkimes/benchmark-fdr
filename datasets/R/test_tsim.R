## ##############################################################################
## test all settings
## ##############################################################################
source("du_tsim.R")
source("funcs_pi0.R")
source("funcs_tstat.R")

## 1. pi0 fixed, icovariate = TRUE
dat <- du_tsim(m = 100, pi0 = 0.5, effect_size = 2, n_sample = 5,
               n_groups = 2, icovariate = TRUE)

## 2. pi0 fixed, icovariate = FALSE
dat <- du_tsim(m = 100, pi0 = 0.5, effect_size = 2, n_sample = 5,
               n_groups = 2, icovariate = FALSE)

## 3. pi0 fixed, icovariate = (function: int -> reals)
dat <- du_tsim(m = 100, pi0 = 0.5, effect_size = 2, n_sample = 5,
               n_groups = 2, icovariate = runif)

## 4. pi0 = (function: int -> reals), icovariate = (function: int -> reals)
dat <- du_tsim(m = 100, pi0 = pi0_smooth1, effect_size = 2, n_sample = 5,
               n_groups = 2, icovariate = runif)

## 5. pi0 = (function: int -> reals), icovariate = logical (ERROR!!) 
dat <- du_tsim(m = 100, pi0 = pi0_smooth1, effect_size = 2, n_sample = 5,
               n_groups = 2, icovariate = TRUE)

## 6. n_sample length 2
dat <- du_tsim(m = 100, pi0 = pi0_smooth1, effect_size = 2, n_sample = c(10, 5),
               n_groups = 2, icovariate = runif)

## 7. effect_size = (function: int -> reals)
dat <- du_tsim(m = 1e5,
               pi0 = .5,
               effect_size = sampler_bimodal,
               icovariate = runif)

