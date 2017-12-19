#!/bin/bash
#SBATCH -J knitChIPseq                 # A single job name for the array
#SBATCH -n 8                        # Number of cores
#SBATCH -N 1                        # All cores on one machine
#SBATCH -p irizarry,serial_requeue  # Partition
#SBATCH --mem 60000                 # Memory request 
#SBATCH -t 0-5:00                   # Maximum execution time (D-HH:MM)
#SBATCH -o render-%j.out  # Standard output
#SBATCH -e render-%j.err  # Standard error
 
module load rstudio
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted. ncores is passed in through environment
# variable SLURM_NTASKS
R -e "rmarkdown::render('~/workspace/cowork/fdr/benchmark-fdr/datasets/ChIP-seq/ChIP-seq.Rmd')"