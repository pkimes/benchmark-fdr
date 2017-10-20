## Submit simulation jobs to slurm cluster
## (tsim: t-test simulation)
##
## Patrick Kimes

M <- 100
ncores <- 10

params <- list(c("pi0", "uniform"),
               c("pi0", "se"),
               c("pi0", "bl"),
               c("n", "uniform"),
               c("n", "se"),
               c("n", "bl"),
               c("esize_fixed", "uniform"),
               c("esize_fixed", "se"),
               c("esize_fixed", "bl"),
               c("esize_random", "uniform"),
               c("esize_random", "se"),
               c("esize_random", "bl"))


for (ip in params) {
    cmd <- paste("sbatch",
                 "-p", "irizarry",
                 "-N", "1",
                 "-n", ncores + 1,
                 "--mem", "64G",
                 "-t", "4-00:00",
                 "-o", "_slurmio/slurm.%x.%j-%a.out",
                 "-e", "_slurmio/slurm.%x.%j-%a.err",
                 "--wrap=\"Rscript tsim-core.R", M, ncores, ip[1], ip[2], "\"")
    system(cmd)
}

