#!/bin/bash
#SBATCH -J knit-sc
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p shared
#SBATCH --mem 125G
#SBATCH -t 0-4:00
#SBATCH -o render-%j.out
#SBATCH -e render-%j.err

module load pandoc
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/pandoc/2.0.2-fasrc01/bin/pandoc/"

# change filename to Rmd to be knitted
# make sure number of cores used in Rmd matches sbatch param -n
R -e "rmarkdown::render('scrnaseq-human-scdd.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-mouse-scdd.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-human-mast.Rmd', clean=FALSE)"
R -e "rmarkdown::render('scrnaseq-mouse-mast.Rmd', clean=FALSE)"



