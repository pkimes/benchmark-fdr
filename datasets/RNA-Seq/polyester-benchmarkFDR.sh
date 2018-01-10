#!/bin/bash

#SBATCH -n 4 # Number of cores
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH -p irizarry,serial_requeue
#SBATCH -e polyester-benchmarkFDR.err
#SBATCH -o polyester-benchmarkFDR.out
#SBATCH --mem=20GB # Memory pool for all cores (see also --mem-per-cpu)

module load centos6/pandoc-1.12.3.3
cd /n/irizarryfs01_backed_up/shicks/projects/benchmark-fdr
Rscript -e 'library(knitr); rmarkdown::render("datasets/RNA-Seq/polyester-benchmarkFDR.Rmd", "html_document")' 
