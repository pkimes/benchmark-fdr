#!/bin/bash
#SBATCH -J knitGWAS                 # A single job name for the array
#SBATCH -n 8                        # Number of cores
#SBATCH -N 1                        # All cores on one machine
#SBATCH -p irizarry,serial_requeue  # Partition
#SBATCH --mem 60000                 # Memory request 
#SBATCH -t 0-2:00                   # Maximum execution time (D-HH:MM)
#SBATCH -o SLURM/render-%j.out  # Standard output
#SBATCH -e SLURM/render-%j.err  # Standard error
 
module load plink
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted. 
# make sure ncores in Rmd matches -n batch param above
R -e "rmarkdown::render('GWAS.Rmd', clean = FALSE)"
