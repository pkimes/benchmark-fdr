#!/bin/bash
#SBATCH -J knitChIPseq
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p irizarry,serial_requeue
#SBATCH --mem 64G
#SBATCH -t 1-05:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err
 
export SLURM_NTASKS
export RSTUDIO_PANDOC="/n/helmod/apps/centos7/Core/rstudio/1.1.453-fasrc01/bin/pandoc"
module load picard

# change filename to Rmd to be knitted. ncores is passed in through environment
# variable SLURM_NTASKS
R -e "rmarkdown::render('ChIPseq-h3k4me3-promoters.Rmd')"
#R -e "rmarkdown::render('ChIPseq-h3k4me3-csaw.Rmd')"
#R -e "rmarkdown::render('ChIPseq-CBP-csaw.Rmd')"
