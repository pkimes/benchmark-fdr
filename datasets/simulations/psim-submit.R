## Submit simulation jobs to slurm cluster
## (psim: direct parameter simulation)
##
## Patrick Kimes

M <- 100
ncores <- 10

params <- list(c("esize_fixed", "uniform", "0.90"),
               c("esize_fixed", "uniform", "0.95"),
               c("esize_fixed", "uniform", "0.99"),
               
               c("esize_fixed", "bl-step-90"),
               c("esize_fixed", "bl-step-95"),
               c("esize_fixed", "bl-cubic-90"),
               c("esize_fixed", "bl-cubic-95"))

               ## c("esize_random_ua", "bl-step-90"),
               ## c("esize_random_ua", "bl-step-95"),
               ## c("esize_random_ua", "bl-cubic-90"),
               ## c("esize_random_ua", "bl-cubic-95"), 

               ## c("esize_random_shift", "bl-step-90"),
               ## c("esize_random_shift", "bl-step-95"),
               ## c("esize_random_shift", "bl-cubic-90"),
               ## c("esize_random_shift", "bl-cubic-95"),

               ## c("altnoise", "bl-step-90"),
               ## c("altnoise", "bl-step-95"),
               ## c("altnoise", "bl-cubic-90"),
               ## c("altnoise", "bl-cubic-95"),

               ## c("pi0", "uniform"),
               ## c("allnull", "uniform"),
               ## c("esize_fixed", "uniform"),
               ## c("esize_fixed", "bl"),
               ## c("esize_fixed", "bl-step-less"),
               ## c("esize_fixed", "bl-step-more"),
               ## c("esize_fixed", "bl-cubic"),
               ## c("esize_random_ua", "uniform"),
               ## c("esize_random_ua", "bl"),
               ## c("esize_random_ua", "bl-step-less"),
               ## c("esize_random_ua", "bl-step-more"),
               ## c("esize_random_ua", "bl-cubic"),
               ## c("esize_random_shift", "uniform"),
               ## c("esize_random_shift", "bl"),
               ## c("esize_random_shift", "bl-step-less"),
               ## c("esize_random_shift", "bl-step-more"),
               ## c("esize_random_shift", "bl-cubic"),
               ## c("altnoise", "uniform"),
               ## c("altnoise", "bl"),
               ## c("altnoise", "bl-step-less"),
               ## c("altnoise", "bl-step-more"),
               ## c("altnoise", "bl-cubic"))

for (ip in params) {
    cmd <- paste("sbatch",
                 "-p", "shared",
                 "-N", "1",
                 "-n", ncores + 1,
                 "--mem", "64G",
                 "-t", "4-00:00",
                 "-o", "_slurmio/slurm.%x.%j-%a.out",
                 "-e", "_slurmio/slurm.%x.%j-%a.err",
                 "--wrap=\"Rscript psim-core.R", M, ncores,
                 paste(ip, collapse = ' '), "\"")
    system(cmd)
}

