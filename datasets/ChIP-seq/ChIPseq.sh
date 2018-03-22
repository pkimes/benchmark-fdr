#!/bin/bash
#SBATCH -J knitChIPseq
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p irizarry,serial_requeue
#SBATCH --mem 16G
#SBATCH -t 0-05:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted. ncores is passed in through environment
# variable SLURM_NTASKS
R -e "rmarkdown::render('ChIPseq-h3k4me3-promoters.Rmd')"
R -e "rmarkdown::render('ChIPseq-h3k4me3-csaw.Rmd')"
