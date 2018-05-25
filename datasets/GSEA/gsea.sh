#!/bin/bash
#SBATCH -J knitGSEA
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p shared,irizarry
#SBATCH --mem 32G
#SBATCH -t 0-6:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
export RSTUDIO_PANDOC="/n/helmod/apps/centos7/Core/rstudio/1.1.453-fasrc01/bin/pandoc"

R -e "rmarkdown::render('gsea-human.Rmd')"
R -e "rmarkdown::render('gsea-mouse.Rmd')"
