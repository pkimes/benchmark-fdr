#$ -cwd
#$ -pe local 1
#$ -l mem_free=2G,h_vmem=2G 

cd /users/shicks1/projects/benchmark-fdr/
Rscript -e 'library(knitr); rmarkdown::render("datasets/RNA-Seq/polyester-replicates-benchmarkFDR.Rmd", "html_document")' 

