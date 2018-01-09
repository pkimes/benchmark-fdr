## Submit simulation jobs to slurm cluster
## (tsim: t-test simulation)
##
## Patrick Kimes

M <- 100
ncores <- 10

params <- list(c("n", "uniform", "0.90"),
               c("n", "uniform", "0.95"),
               c("n", "uniform", "0.99"),
               c("n", "se", "0.90"),
               c("n", "se", "0.95"),
               c("n", "se", "0.99"),
               c("n", "bl-step-90"),
               c("n", "bl-step-95"),
               c("n", "bl-cubic-90"),
               c("n", "bl-cubic-95"),

               c("esize_fixed", "uniform", "0.90"),
               c("esize_fixed", "uniform", "0.95"),
               c("esize_fixed", "uniform", "0.99"),
               c("esize_fixed", "se", "0.90"),
               c("esize_fixed", "se", "0.95"),
               c("esize_fixed", "se", "0.99"),
               c("esize_fixed", "bl-step-90"),
               c("esize_fixed", "bl-step-95"),
               c("esize_fixed", "bl-cubic-90"),
               c("esize_fixed", "bl-cubic-95"),

               c("esize_random", "uniform", "0.90"),
               c("esize_random", "uniform", "0.95"),
               c("esize_random", "uniform", "0.99"),
               c("esize_random", "se", "0.90"),
               c("esize_random", "se", "0.95"),
               c("esize_random", "se", "0.99"),
               c("esize_random", "bl-step-90"),
               c("esize_random", "bl-step-95"),
               c("esize_random", "bl-cubic-90"),
               c("esize_random", "bl-cubic-95"),
               
               c("pi0", "uniform"),
               c("pi0", "se"),
               c("n", "uniform"),
               c("n", "se"),
               c("n", "bl"),
               c("n", "bl-step-less"),
               c("n", "bl-step-more"),
               c("n", "bl-cubic"),
               c("esize_fixed", "uniform"),
               c("esize_fixed", "se"),
               c("esize_fixed", "bl"),
               c("esize_fixed", "bl-step-less"),
               c("esize_fixed", "bl-step-more"),
               c("esize_fixed", "bl-cubic"),
               c("esize_random", "uniform"),
               c("esize_random", "se"),
               c("esize_random", "bl"),
               c("esize_random", "bl-step-less"),
               c("esize_random", "bl-step-more"),
               c("esize_random", "bl-cubic"))


for (ip in params) {
    cmd <- paste("sbatch",
                 "-p", "shared",
                 "-N", "1",
                 "-n", ncores + 1,
                 "--mem", "64G",
                 "-t", "4-00:00",
                 "-o", "_slurmio/slurm.%x.%j-%a.out",
                 "-e", "_slurmio/slurm.%x.%j-%a.err",
                 "--wrap=\"Rscript tsim-core.R", M, ncores,
                 paste(ip, collapse = ' '), "\"")
    system(cmd)
}

