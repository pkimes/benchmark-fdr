# Submit simulation jobs to slurm cluster
#
# Patrick Kimes

M <- 100
ncores <- 10

for (sidx in 1:2) {
    cmd <- paste("sbatch",
                 "-p", "irizarry",
                 "-N", "1",
                 "-n", ncores + 1,
                 "--mem", "64G",
                 "-t", "4-00:00",
                 "-o", "_slurmio/slurm.%x.%j-%a.out",
                 "-e", "_slurmio/slurm.%x.%j-%a.err",
                 "--wrap=\"Rscript sim-core.R", M, ncores, sidx, "\"")
    system(cmd)
}

