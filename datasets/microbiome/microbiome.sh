#!/bin/bash
#SBATCH -J knit16s
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p shared,irizarry
#SBATCH --mem 32G
#SBATCH -t 0-8:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
##module load rstudio
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"

R -e "rmarkdown::render('microbiome-baxter.Rmd')"
R -e "rmarkdown::render('microbiome-goodrich.Rmd')"
R -e "rmarkdown::render('microbiome-schubert.Rmd')"
R -e "rmarkdown::render('microbiome-enigma.Rmd')"
R -e "rmarkdown::render('microbiome-papa.Rmd')"
