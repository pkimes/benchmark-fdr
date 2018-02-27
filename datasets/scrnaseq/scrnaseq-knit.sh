#!/bin/bash
#SBATCH -J knit-sc4                  # A single job name for the array
#SBATCH -n 6                         # Number of cores
#SBATCH -N 1                         # All cores on one machine
#SBATCH -p irizarry,serial_requeue   # Partition
#SBATCH --mem 200G                   # Memory request 
#SBATCH -t 0-20:00                   # Maximum execution time (D-HH:MM)
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
module load rstudio
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted. ncores is passed in through environment
# variable SLURM_NTASKS
R -e "rmarkdown::render('scrnaseq-human-scdd.Rmd')"
R -e "rmarkdown::render('scrnaseq-human-mast.Rmd')"
R -e "rmarkdown::render('scrnaseq-mouse-scdd.Rmd')"
R -e "rmarkdown::render('scrnaseq-mouse-mast.Rmd')"
