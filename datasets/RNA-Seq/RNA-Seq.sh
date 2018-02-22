#!/bin/bash
#SBATCH -J knitRNAseq
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p shared
#SBATCH --mem 32G
#SBATCH -t 0-6:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
##module load rstudio
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted. ncores is passed in through environment
# variable SLURM_NTASKS
R -e "rmarkdown::render('RNA-Seq.Rmd')"
