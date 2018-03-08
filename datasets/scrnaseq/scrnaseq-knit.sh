#!/bin/bash
#SBATCH -J knit-sc                  # A single job name for the array
#SBATCH -n 8                        # Number of cores
#SBATCH -N 1                        # All cores on one machine
#SBATCH -p shared                   # Partition
#SBATCH --mem 125G                  # Memory request 
#SBATCH -t 0-6:00                   # Maximum execution time (D-HH:MM)
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted
# make sure number of cores used in Rmd matches sbatch param -n
R -e "rmarkdown::render('scrnaseq-mouse-scdd.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-mouse-mast.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-human-scdd.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-human-mast.Rmd', clean=FALSE)"
