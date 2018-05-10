#!/bin/bash
#SBATCH -J knitGSEA
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p shared,irizarry
#SBATCH --mem 32G
#SBATCH -t 0-6:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
##module load rstudio
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

R -e "rmarkdown::render('gsea-human.Rmd')"
R -e "rmarkdown::render('gsea-mouse.Rmd')"
