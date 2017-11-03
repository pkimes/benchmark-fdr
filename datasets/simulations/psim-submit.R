## Submit simulation jobs to slurm cluster
## (psim: direct parameter simulation)
##
## Patrick Kimes

M <- 100
ncores <- 10

params <- list(c("pi0", "uniform"),
               c("esize_fixed", "uniform"),
               c("esize_fixed", "bl"),
               c("esize_random_ua", "uniform"),
               c("esize_random_ua", "bl"),
               c("esize_random_shift", "uniform"),
               c("esize_random_shift", "bl"),
               c("altnoise", "uniform"),
               c("altnoise", "bl"))

for (ip in params) {
    cmd <- paste("sbatch",
                 "-p", "irizarry",
                 "-N", "1",
                 "-n", ncores + 1,
                 "--mem", "64G",
                 "-t", "4-00:00",
                 "-o", "_slurmio/slurm.%x.%j-%a.out",
                 "-e", "_slurmio/slurm.%x.%j-%a.err",
                 "--wrap=\"Rscript psim-core.R", M, ncores, ip[1], ip[2], "\"")
    system(cmd)
}

