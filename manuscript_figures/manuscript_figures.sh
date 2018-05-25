#!/bin/bash
#SBATCH -J figs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p irizarry
#SBATCH --mem 50G
#SBATCH -t 0-1:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err

export RSTUDIO_PANDOC="/n/helmod/apps/centos7/Core/rstudio/1.1.453-fasrc01/bin/pandoc"

# change filename to Rmd to be knitted
# make sure number of cores used in Rmd matches sbatch param -n
R -e "rmarkdown::render('manuscript_figures.Rmd', clean=TRUE)"


